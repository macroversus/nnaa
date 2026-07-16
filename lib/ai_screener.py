from __future__ import annotations

import json
import logging
import os
import re
import time
from typing import Any, Dict, List, Optional, Set, Tuple

import requests

from deduplicator import normalize_doi, sanitize_csv_text
from nnaa_queries import ALL_NNAA_TRACKS

logger = logging.getLogger(__name__)

AI_EXTRA_COLUMNS = ("ai_nnaa_category", "ai_synthesis_method", "ai_has_pathway_or_synthesis", "ai_compound")

_AI_BASE_COLUMNS = (
    "ai_article_type",
    "ai_domain_relevant",
    "ai_has_experiment",
    "ai_pass",
    "ai_rationale_zh",
    "lit_type_zh",
)

_TRACK_ORDER = {track: idx for idx, track in enumerate(ALL_NNAA_TRACKS)}
_TRACK_ORDER[""] = len(ALL_NNAA_TRACKS)

_VALID_NNAA_CATEGORIES = frozenset(
    {"pathway", "enzymatic", "fermentation", "chemical", "hybrid", "gce", "both", "unclear", "neither"}
)
_VALID_SYNTHESIS_METHODS = frozenset(
    {"enzymatic", "fermentation", "chemical", "hybrid", "multiple", "unclear", "none"}
)


_LIT_TYPE_ZH_MAP = {
    "review": "综述",
    "research_article": "研究论文",
    "editorial_comment_letter": "评论/通讯",
    "other": "其他",
}


def _normalize_nnaa_category(val):
    s = str(val or "").strip().lower()
    if s in _VALID_NNAA_CATEGORIES:
        return s
    return ""


def _normalize_synthesis_method(val):
    s = str(val or "").strip().lower()
    if s in _VALID_SYNTHESIS_METHODS:
        return s
    return ""


def _normalize_yes_no_unclear(val):
    s = str(val or "").strip().lower()
    if s in ("yes", "no", "unclear"):
        return s
    return ""


def _set_ai_field_defaults(r):
    for k in _AI_BASE_COLUMNS:
        r.setdefault(k, "")
    for k in AI_EXTRA_COLUMNS:
        r.setdefault(k, "")


def _apply_ai_judgment_fields(r, j):
    r["ai_nnaa_category"] = _normalize_nnaa_category(j.get("nnaa_category"))
    r["ai_synthesis_method"] = _normalize_synthesis_method(j.get("synthesis_method"))
    r["ai_has_pathway_or_synthesis"] = _normalize_yes_no_unclear(j.get("has_pathway_or_synthesis"))
    r["ai_compound"] = sanitize_csv_text(j.get("compound_name") or "")[:200]
    # 生成中文文献类型标签
    art_type = str(r.get("ai_article_type") or "").strip().lower()
    r["lit_type_zh"] = _LIT_TYPE_ZH_MAP.get(art_type, "其他")


def _strip_json_fence(text: str) -> str:
    t = text.strip()
    if t.startswith("```"):
        t = re.sub(r"^```(?:json)?\s*", "", t, flags=re.IGNORECASE)
        t = re.sub(r"\s*```$", "", t)
    return t.strip()


def _call_chat(
    api_base: str,
    api_key: str,
    model: str,
    system: str,
    user: str,
    timeout: int = 120,
) -> str:
    base = api_base.rstrip("/")
    url = f"{base}/chat/completions"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    body: Dict[str, Any] = {
        "model": model,
        "temperature": 0.2,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
    }
    resp = requests.post(url, headers=headers, json=body, timeout=timeout)
    resp.raise_for_status()
    data = resp.json()
    choices = data.get("choices") or []
    if not choices:
        raise ValueError("API 返回无 choices")
    msg = choices[0].get("message") or {}
    return (msg.get("content") or "").strip()


def _as_str(val: Any, max_len: int = 0) -> str:
    s = sanitize_csv_text(val)
    return s[:max_len] if max_len else s


def _build_user_payload(batch: List[Dict]) -> str:
    slim = []
    for rec in batch:
        slim.append(
            {
                "doi": _as_str(rec.get("doi")),
                "title": _as_str(rec.get("title"), 800),
                "journal": _as_str(rec.get("journal")),
                "pub_date": _as_str(rec.get("pub_date")),
                "source": _as_str(rec.get("source")),
                "nnaa_track": _as_str(rec.get("nnaa_track") or rec.get("enzyme_track")),
                "pmid": _as_str(rec.get("pmid")),
                "publication_types": _as_str(rec.get("publication_types"), 500),
                "abstract": _as_str(rec.get("abstract"), 3500),
            }
        )
    return json.dumps(slim, ensure_ascii=False)


def _batch_dois(batch: List[Dict]) -> Set[str]:
    return {normalize_doi(r.get("doi") or "") for r in batch if normalize_doi(r.get("doi") or "")}


def _record_track(rec: Dict) -> str:
    track = str(rec.get("nnaa_track") or rec.get("enzyme_track") or "").strip().lower()
    if track in ALL_NNAA_TRACKS:
        return track
    return ""


def _sort_for_screening(records: List[Dict]) -> List[Dict]:
    return sorted(
        records,
        key=lambda r: (
            _TRACK_ORDER.get(_record_track(r), 9),
            normalize_doi(r.get("doi") or ""),
        ),
    )


_BASE_SYSTEM = """你是非天然氨基酸（NNAA / unnatural amino acid / noncanonical amino acid）领域文献筛选助手。
输入为 JSON 数组，每项有一条文献的 doi、title、journal、pub_date、source、nnaa_track、pmid、publication_types、abstract。
请对**每一条**输出一个判断对象，组成 JSON 数组，顺序与输入一致，且每个对象必须包含字段 doi（与输入一致）。

══════════════════════════════════════════════
【第一步：化合物识别——它是氨基酸吗？】
══════════════════════════════════════════════
首先判断文献的核心研究化合物是否属于氨基酸类（α-氨基酸、β-氨基酸、N-取代氨基酸、氨基酸衍生物）。
如果核心研究对象不是氨基酸，直接 domain_relevant=false，不需要继续判断。

明确**不是**氨基酸的化合物（这些出现时直接排除）：
  ✗ 糖类：葡萄糖、蔗糖、淀粉、壳聚糖、透明质酸、麦芽糖、纤维素等
  ✗ 萜类：紫杉醇、青蒿素、倍半萜、单萜、二萜、三萜等
  ✗ 生物碱：长春碱、喜树碱、可卡因等（虽含氮但非氨基酸）
  ✗ 有机酸/聚酮：衣康酸、戊二酸、对香豆酸、丙烯酸、琥珀酸、己二酸等
  ✗ 维生素：维生素B12、维生素C、叶酸等
  ✗ 核苷/核苷酸：ATP、NAD、CoA 等
  ✗ 脂肪酸/聚羟基脂肪酸（PHA/PHB）：癸二酸、油酸等
  ✗ 醇类：1,3-丙二醇、丁二醇等
  ✗ 芳香族化合物（非氨基酸骨架）：苯乙烯、肉桂酸、对苯二酚等
  注意：以上化合物在代谢工程论文中频繁出现，但不属于 NNAA 范畴。

══════════════════════════════════════════════
【第二步：是天然蛋白质氨基酸还是非天然氨基酸？】
══════════════════════════════════════════════
20 种标准蛋白质氨基酸：Gly, Ala, Val, Leu, Ile, Pro, Phe, Trp, Met, Ser, Thr, Cys, Tyr, His, Asp, Glu, Asn, Gln, Lys, Arg。
以下情形即使研究的是天然氨基酸本身，也必须 domain_relevant=false：
  ✗ 谷氨酸（glutamic acid / monosodium glutamate）的发酵生产
  ✗ 赖氨酸（L-lysine）的工业发酵（作为饲料添加剂）
  ✗ 色氨酸、苏氨酸等天然氨基酸的普通代谢研究
  ✗ 仅研究天然氨基酸的代谢通路作为背景信息

例外：下列"天然"氨基酸在 NNAA 研究中有特殊地位，若文献以其为 NNAA 工具/底物则相关：
  △ 鸟氨酸 (ornithine)、瓜氨酸 (citrulline) — 作为 NNAA 代谢前体/中间体
  △ 硒代半胱氨酸 (selenocysteine)、吡咯赖氨酸 (pyrrolysine) — 第21/22种天然 AA，GCE 研究核心
  △ 羟脯氨酸 (hydroxyproline)、磷酸丝氨酸 (phosphoserine) — 翻译后修饰研究

══════════════════════════════════════════════
【第三步：NNAA 研究范畴认定（需通过前两步）】
══════════════════════════════════════════════
认定为相关（domain_relevant=true）的情形：
  ✓ NNAA 的化学合成、不对称合成、立体选择性制备
  ✓ NNAA 的酶法合成（转氨酶、氨酰-tRNA 合成酶催化等）
  ✓ NNAA 的发酵/代谢工程生产（细胞工厂）
  ✓ NNAA 的生物合成/代谢通路研究（途径重构、前体供给等）
  ✓ 遗传密码扩展（GCE）：aaRS/tRNA 系统介导的 NNAA 蛋白质整合
  ✓ NNAA 在多肽、拟肽、肽类药物中的应用（以 NNAA 为核心）
  ✓ NNAA 的结构、性质、生物活性研究
  ✓ D-氨基酸的合成、代谢、转化

典型 NNAA 化合物（见到即相关）：
  norvaline, norleucine, tert-leucine, Aib, alpha-aminoisobutyric acid,
  cyclohexylalanine, 4-fluorophenylalanine, naphthylalanine, biphenylalanine,
  4-hydroxyproline, pipecolic acid, citrulline (as NNAA), homoarginine,
  5-hydroxytryptophan, D-amino acids (D-Ala, D-Phe etc.), beta-amino acids,
  p-azidophenylalanine, p-benzoylphenylalanine, p-acetylphenylalanine,
  BocK, AllocK, propargyllysine, phosphoserine, 3-nitrotyrosine, pyrrolysine,
  azidonorleucine, DOPA (3,4-dihydroxyphenylalanine)

══════════════════════════════════════════════
【输出字段（JSON 数组，每项必含 doi）】
══════════════════════════════════════════════
1) article_type: "review" | "research_article" | "editorial_comment_letter" | "other"

2) domain_relevant: true/false — 按以上三步判断，NNAA 是否为核心研究对象。

3) has_experiment: "yes" | "no" | "unclear"

4) has_pathway_or_synthesis: "yes" | "no" | "unclear"
   文献是否明确包含以下任一可供数据库使用的内容：
   → yes：有具体 NNAA 代谢/生物合成通路（酶、基因、底物、产物）
   → yes：有具体 NNAA 化学合成路线（反应步骤、试剂、收率、ee值）
   → yes：有 NNAA 酶法制备工艺（酶种类、反应条件、转化率）
   → yes：有 NNAA 发酵/代谢工程生产方案（菌株改造、滴度、产率）
   → yes：有 NNAA 前体/中间体的合成分析
   → unclear：摘要提示有上述内容但细节不足以判断
   → no：文献只是"使用"已知 NNAA（蛋白标记/click chemistry/生物成像等），未报告如何制备
   → no：只研究 NNAA 结构、性质、生物活性，无合成/通路信息

5) pass_filter: true/false
   按文章类型分两套规则：

   【综述（review）】满足以下两点即为 true：
   ① domain_relevant = true
   ② has_pathway_or_synthesis = yes 或 unclear（综述涵盖 NNAA 通路/合成方法，对数据库构建有参考价值）

   【研究论文（research_article）】同时满足以下三点才为 true：
   ① domain_relevant = true
   ② has_experiment = yes 或 unclear
   ③ has_pathway_or_synthesis = yes 或 unclear

   【评论/通讯/其他（editorial_comment_letter / other）】一律为 false

   注意：仅"使用" NNAA（蛋白标记/click反应/成像等）而无制备/通路内容 → pass_filter 必须为 false，无论哪种文章类型

6) rationale_zh: 1-2句中文理由。pass=false 时须明确说明原因：
   若非 NNAA 化合物 → 说明实际研究的是什么
   若纯应用 → 说明"仅使用已知ncAA做X，未报告合成或通路"
   若缺乏通路/合成内容 → 说明文献内容是什么

7) nnaa_category: "pathway"|"enzymatic"|"fermentation"|"chemical"|"hybrid"|"gce"|"both"|"unclear"|"neither"

8) synthesis_method: "enzymatic"|"fermentation"|"chemical"|"hybrid"|"multiple"|"unclear"|"none"

9) compound_name: 文献核心 NNAA 化合物名（英文通用名或缩写，逗号分隔，最多5个；非 NNAA 填 ""）

只输出 JSON 数组，不要 Markdown 代码围栏，不要其它文字。"""

_TRACK_FOCUS = {
    "pathway": """
当前轨道：非天然氨基酸代谢通路（pathway）。
重点：NNAA 生物合成/降解通路、NNAA 代谢工程、途径重构、前体供给、转运、调控、合成生物学底盘中的 NNAA 代谢网络。
严格过滤：文献主题必须是 NNAA 的通路。以下情形绝对排除：
  - 研究天然氨基酸（谷氨酸、赖氨酸等）代谢通路（非作为 NNAA 前体）
  - 研究有机酸、脂肪酸、萜类、维生素等非氨基酸化合物的代谢通路
  - 通路只作为背景信息，核心研究的化合物不是 NNAA""",
    "enzymatic": """
当前轨道：酶法合成（enzymatic）。
重点：转氨酶/氨转移酶、氨酰-tRNA 合成酶、酶级联、体外翻译、酶促偶联等 NNAA 酶法制备。
严格过滤：酶的底物或产物必须是 NNAA；仅研究酶的结构/机理但不制备 NNAA 则排除。""",
    "fermentation": """
当前轨道：生物发酵（fermentation）。
重点：大肠杆菌/酵母/谷氨酸棒杆菌等微生物细胞工厂、发酵代谢流改造、全细胞催化生产 NNAA。
严格过滤：发酵目标产物必须是 NNAA；生产天然氨基酸（谷氨酸、赖氨酸等）的发酵工艺绝对排除。""",
    "chemical": """
当前轨道：化学合成（chemical）。
重点：NNAA 不对称合成、全合成、有机催化、光/redox 驱动、保护基/手性诱导等化学路线。
严格过滤：合成目标必须是 NNAA；合成其他手性化合物（糖、萜类等）不相关。""",
    "hybrid": """
当前轨道：酶-化学联用（hybrid）。
重点：chemoenzymatic、one-pot 多步、酶促步骤与化学步骤串联/并联的 NNAA 制备策略。""",
    "gce": """
当前轨道：遗传密码扩展（GCE）。
重点：aaRS/tRNA 系统介导的 NNAA 定点整合到蛋白质、amber suppression（TAG/TGA 密码子抑制）、
       遗传密码扩展的 ncAA（p-azidophenylalanine、p-benzoylphenylalanine、propargyllysine 等）、
       aaRS 工程化（MjTyrRS、PylRS、EcTyrRS 等）、正交 tRNA/aaRS 对、GCE 应用（蛋白质标记、交联、生物正交化学）。
注意：GCE 论文的核心化合物是 ncAA，即使文章重点在蛋白质工程，也应认定为 domain_relevant=true。""",
}

_DEFAULT_FOCUS = """
若 nnaa_track 缺失，综合判断 NNAA 代谢通路与各合成路线（酶法/发酵/化学/联用）相关主题。"""


def build_system_prompt(ai_cfg: dict, nnaa_track: str) -> str:
    track = (nnaa_track or "").strip().lower()
    custom_key = f"system_prompt_{track}" if track in ALL_NNAA_TRACKS else ""
    custom = (ai_cfg.get(custom_key) or "").strip() if custom_key else ""
    if custom:
        system = custom
    else:
        system = _BASE_SYSTEM
        if track in _TRACK_FOCUS:
            system += _TRACK_FOCUS[track]
        else:
            system += _DEFAULT_FOCUS

    extra_key = f"extra_instructions_{track}" if track in ALL_NNAA_TRACKS else ""
    extra = (ai_cfg.get(extra_key) or ai_cfg.get("extra_instructions") or "").strip()
    if extra:
        system += "\n\n补充说明：\n" + extra
    return system


def screen_records(
    records: List[Dict],
    ai_cfg: dict,
    *,
    state_dir: Optional[str] = None,
    resume: bool = False,
    clear_checkpoint: bool = False,
) -> Tuple[List[Dict], Dict[str, int]]:
    stats = {
        "ai_input": len(records),
        "ai_batches": 0,
        "ai_batches_skipped": 0,
        "ai_dropped": 0,
        "ai_kept": 0,
        "ai_errors": 0,
    }
    if not records or not ai_cfg.get("enabled", False):
        out: List[Dict] = []
        for rec in records:
            r = dict(rec)
            _set_ai_field_defaults(r)
            out.append(r)
        return out, stats

    api_key_env = ai_cfg.get("api_key_env", "NNAA_FETCH_OPENAI_API_KEY")
    api_key = os.environ.get(api_key_env, "").strip() or (ai_cfg.get("api_key") or "").strip()
    api_base = (
        os.environ.get("NNAA_FETCH_OPENAI_BASE", "").strip()
        or ai_cfg.get("api_base", "https://api.openai.com/v1")
    ).rstrip("/")
    model = ai_cfg.get("model", "gpt-4o-mini")
    batch_size = max(1, int(ai_cfg.get("batch_size", 6)))
    only_keep = ai_cfg.get("only_keep_passed", True)
    keep_on_batch_error = ai_cfg.get("keep_on_batch_error", True)
    interval = float(ai_cfg.get("request_interval_sec", 0.4))
    fail_missing = ai_cfg.get("fail_if_missing_key", False)
    checkpoint_enabled = bool(state_dir and ai_cfg.get("checkpoint_enabled", True))

    if not api_key:
        msg = f"AI 筛选已启用但未设置 API Key（环境变量 {api_key_env} 或配置 api_key）"
        if fail_missing:
            raise RuntimeError(msg)
        logger.warning(f"{msg}，跳过 AI，全部保留。")
        out = []
        for rec in records:
            r = dict(rec)
            r["ai_article_type"] = ""
            r["ai_domain_relevant"] = ""
            r["ai_has_experiment"] = ""
            r["ai_pass"] = ""
            r["ai_rationale_zh"] = "skipped_no_api_key"
            for k in AI_EXTRA_COLUMNS:
                r[k] = ""
            out.append(r)
        stats["ai_kept"] = len(out)
        return out, stats

    ordered_records = _sort_for_screening(records)
    signature = ""
    by_doi: Dict[str, Dict] = {}
    batch_failed: Set[str] = set()

    if checkpoint_enabled:
        from ai_checkpoint import (
            clear_ai_checkpoint,
            load_ai_checkpoint,
            records_signature,
            save_ai_checkpoint,
        )

        signature = records_signature(ordered_records)
        if clear_checkpoint:
            clear_ai_checkpoint(state_dir)
        elif not resume:
            existing = load_ai_checkpoint(state_dir)
            if existing and existing.get("signature") != signature:
                clear_ai_checkpoint(state_dir)
        if resume:
            cp = load_ai_checkpoint(state_dir, expected_signature=signature)
            if cp:
                by_doi = dict(cp.get("results_by_doi") or {})
                batch_failed = set(cp.get("batch_failed_dois") or [])

    processed = set(by_doi.keys()) | batch_failed

    for i in range(0, len(ordered_records), batch_size):
        batch = ordered_records[i : i + batch_size]
        stats["ai_batches"] += 1
        batch_dois = _batch_dois(batch)

        if batch_dois and batch_dois.issubset(processed):
            stats["ai_batches_skipped"] += 1
            continue

        track = _record_track(batch[0]) if batch else ""
        system = build_system_prompt(ai_cfg, track)
        user = _build_user_payload(batch)
        try:
            raw = _call_chat(api_base, api_key, model, system, user)
            parsed = json.loads(_strip_json_fence(raw))
            if not isinstance(parsed, list):
                raise ValueError("顶层不是数组")
            for item in parsed:
                if not isinstance(item, dict):
                    continue
                doi = normalize_doi(item.get("doi") or "")
                if not doi:
                    continue
                by_doi[doi] = item
                batch_failed.discard(doi)
        except Exception as e:
            stats["ai_errors"] += 1
            logger.warning(f"AI 批次解析失败 (batch {stats['ai_batches']}, track={track or 'mixed'}): {e}")
            if keep_on_batch_error:
                for rec in batch:
                    d = normalize_doi(rec.get("doi") or "")
                    if d:
                        batch_failed.add(d)
                        by_doi.pop(d, None)
        else:
            for d in batch_dois:
                if d not in by_doi and d not in batch_failed:
                    logger.debug(f"AI 批次 {stats['ai_batches']} 未返回 DOI: {d}")

        processed = set(by_doi.keys()) | batch_failed
        if checkpoint_enabled:
            from ai_checkpoint import save_ai_checkpoint

            save_ai_checkpoint(
                state_dir,
                signature,
                len(ordered_records),
                batch_size,
                by_doi,
                batch_failed,
                stats["ai_batches"],
            )
        time.sleep(interval)

    if checkpoint_enabled:
        from ai_checkpoint import clear_ai_checkpoint

        all_dois = {
            normalize_doi(r.get("doi") or "")
            for r in ordered_records
            if normalize_doi(r.get("doi") or "")
        }
        if all_dois.issubset(processed):
            clear_ai_checkpoint(state_dir)
            logger.info("AI 筛选全部完成，已清除 checkpoint")

    if stats["ai_batches_skipped"]:
        logger.info(
            f"AI checkpoint 跳过已完成批次: {stats['ai_batches_skipped']}/{stats['ai_batches']}"
        )

    merged: List[Dict] = []
    for rec in records:
        r = dict(rec)
        doi = normalize_doi(r.get("doi") or "")
        if doi in batch_failed:
            r["ai_article_type"] = ""
            r["ai_domain_relevant"] = ""
            r["ai_has_experiment"] = ""
            r["ai_pass"] = "true"
            r["ai_rationale_zh"] = "batch_api_error_kept"
            for k in AI_EXTRA_COLUMNS:
                r[k] = ""
            merged.append(r)
            stats["ai_kept"] += 1
            continue

        j = by_doi.get(doi, {})
        r["ai_article_type"] = str(j.get("article_type", "") or "")
        dr = j.get("domain_relevant", "")
        r["ai_domain_relevant"] = str(dr).lower() if not isinstance(dr, bool) else ("true" if dr else "false")
        r["ai_has_experiment"] = str(j.get("has_experiment", "") or "")
        pf = j.get("pass_filter", False)
        r["ai_pass"] = "true" if pf else "false"
        r["ai_rationale_zh"] = (j.get("rationale_zh") or "").replace("\n", " ").strip()[:500]
        if not j:
            r["ai_rationale_zh"] = (r["ai_rationale_zh"] or "ai_missing_judgment").strip()
        _apply_ai_judgment_fields(r, j)
        for k in AI_EXTRA_COLUMNS:
            r.setdefault(k, "")

        if only_keep and r["ai_pass"] != "true":
            stats["ai_dropped"] += 1
            continue
        merged.append(r)
        stats["ai_kept"] += 1

    return merged, stats

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
    {"pathway", "enzymatic", "fermentation", "chemical", "hybrid", "both", "unclear", "neither"}
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

【核心判定原则】
文献的**核心研究对象**必须是非天然氨基酸（NNAA）本身，才能认定为相关。
- NNAA 范围：不存在于标准蛋白质编码体系的氨基酸，包括 norvaline、norleucine、tert-leucine、Aib（alpha-aminoisobutyric acid）、homoarginine、citrulline、ornithine（作为 NNAA 前体时）、4-hydroxyproline、pipecolic acid、D-氨基酸、fluorophenylalanine、naphthylalanine、5-hydroxytryptophan 等，以及以上化合物的合成/代谢/应用。
- 以下情形即使涉及氨基酸、代谢工程、发酵也必须标 domain_relevant=false：
  * 研究对象为天然氨基酸（甘氨酸、丙氨酸、谷氨酸、谷氨酰胺、赖氨酸（本身非NNAA）等）的普通生产
  * 研究对象为脂肪酸、有机酸（谷氨酸、戊二酸、衣康酸、对香豆酸等）、多糖、生物碱、萜类等**非氨基酸**的合成/代谢
  * 研究主题为植物发育、病原体互作、基因组学、转录组学，仅提及氨基酸或代谢通路作为背景
  * 文章虽然研究"手性化合物"或"生物催化"但未明确涉及 NNAA 化合物

判断字段：
1) article_type: "review" | "research_article" | "editorial_comment_letter" | "other"
2) domain_relevant: true/false — 文献核心是否聚焦 NNAA 的合成、代谢通路或在蛋白质/肽链中的应用；不符合上述核心原则则 false。
3) has_experiment: "yes" | "no" | "unclear" — 是否包含实验（体外酶促、发酵、化学合成、细胞/动物实验等）。
4) pass_filter: true/false — 同时满足以下三点才能为 true：① 非 review/editorial；② domain_relevant=true；③ has_experiment 为 yes 或 unclear（且标题/摘要强烈暗示有实验）。
5) rationale_zh: 一两句中文理由，若拒绝请指出文献实际研究的是什么化合物/主题。
6) nnaa_category: "pathway" | "enzymatic" | "fermentation" | "chemical" | "hybrid" | "both" | "unclear" | "neither"
   — pathway=NNAA 生物合成/降解通路；enzymatic=酶法制备 NNAA；fermentation=细胞工厂发酵生产 NNAA；chemical=化学路线合成 NNAA；hybrid=酶-化学联用；both=同时涵盖通路+合成；neither=非 NNAA 主题。
7) synthesis_method: "enzymatic" | "fermentation" | "chemical" | "hybrid" | "multiple" | "unclear" | "none"
   — 文献重点报道的 NNAA 合成/制备方法；纯通路研究无合成可标 none；domain_relevant=false 时填 "none"。
8) has_pathway_or_synthesis: "yes" | "no" | "unclear" — 是否明确报道 NNAA 代谢通路或具体合成路线/制备工艺。
9) compound_name: 文献中明确报道的**非天然**氨基酸化合物名称（用最通用英文名或缩写，多个用英文逗号分隔，最多5个）。
   若文献研究的不是 NNAA，或未明确涉及具体 NNAA 化合物则填 ""。
   示例："norvaline, 4-fluorophenylalanine" 或 "Aib" 或 "4-hydroxyproline, pipecolic acid"

只输出 JSON 数组，不要 Markdown 代码围栏，不要其它文字。"""

_TRACK_FOCUS = {
    "pathway": """
当前轨道：非天然氨基酸代谢通路（pathway）。
重点：NNAA 生物合成/降解通路、NNAA 代谢工程、途径重构、前体供给、转运、调控、合成生物学底盘中的 NNAA 代谢网络。
注意：文献主题必须是 NNAA（非天然氨基酸）的通路，研究天然氨基酸代谢、脂肪酸、有机酸或其它化合物的通路不算相关。""",
    "enzymatic": """
当前轨道：酶法合成（enzymatic）。
重点：转氨酶/氨转移酶、氨酰-tRNA 合成酶、酶级联、体外翻译、酶促偶联等 NNAA 酶法制备。""",
    "fermentation": """
当前轨道：生物发酵（fermentation）。
重点：大肠杆菌/酵母/谷氨酸棒杆菌等微生物细胞工厂、发酵代谢流改造、全细胞催化生产 NNAA。""",
    "chemical": """
当前轨道：化学合成（chemical）。
重点：NNAA 不对称合成、全合成、有机催化、光/redox 驱动、保护基/手性诱导等化学路线。""",
    "hybrid": """
当前轨道：酶-化学联用（hybrid）。
重点：chemoenzymatic、one-pot 多步、酶促步骤与化学步骤串联/并联的 NNAA 制备策略。""",
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

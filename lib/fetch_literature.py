#!/usr/bin/env python3

import csv
import os
import sys
import json
import logging
import argparse
from datetime import datetime, date, timedelta
from pathlib import Path
from typing import List, Dict, Optional

import yaml
import pandas as pd

LIB_DIR = Path(__file__).parent.resolve()
ROOT_DIR = LIB_DIR.parent
sys.path.insert(0, str(ROOT_DIR))
sys.path.insert(0, str(LIB_DIR))

from pubmed_fetcher import PubMedFetcher
from crossref_fetcher import CrossRefFetcher
from deduplicator import (
    load_existing_dois_from_config,
    deduplicate_new_records,
    sanitize_record_for_csv,
)
from ai_screener import screen_records
from nnaa_queries import (
    ALL_NNAA_TRACKS,
    resolve_nnaa_tracks,
    resolve_pubmed_queries,
    resolve_crossref_queries,
    resolve_crossref_query_batches,
    merged_crossref_cfg,
)

STATE_FILE = "last_run.json"


def setup_logging(log_dir: str, verbose: bool = False):
    log_dir_path = Path(log_dir)
    log_dir_path.mkdir(parents=True, exist_ok=True)

    log_file = log_dir_path / f"fetch_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, encoding="utf-8"),
        ],
    )
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    return logging.getLogger("main")


def load_state(state_dir: str) -> dict:
    path = Path(state_dir) / STATE_FILE
    if path.exists():
        try:
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def save_state(state_dir: str, state: dict):
    Path(state_dir).mkdir(parents=True, exist_ok=True)
    path = Path(state_dir) / STATE_FILE
    with open(path, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2, ensure_ascii=False)


def validate_config(cfg: dict, logger: logging.Logger) -> List[str]:
    issues: List[str] = []
    paths_cfg = cfg.get("paths") or {}

    master_csv = (paths_cfg.get("master_csv") or "").strip()
    if master_csv and not Path(master_csv).exists():
        issues.append(f"主库 CSV 不存在: {master_csv}")

    for src in paths_cfg.get("extra_doi_csvs") or []:
        p = (src or {}).get("path") or ""
        if p and not Path(p).exists():
            issues.append(f"附加 DOI 文件不存在: {p}")

    pubmed_cfg = cfg.get("pubmed") or {}
    if pubmed_cfg.get("enabled", True):
        has_pm = any(
            resolve_pubmed_queries(cfg, track)
            for track in ALL_NNAA_TRACKS
        ) or bool(pubmed_cfg.get("queries"))
        if has_pm:
            email = (pubmed_cfg.get("email") or "").strip()
            if not email or "example.com" in email.lower():
                issues.append("PubMed email 未配置或为占位符，请在 config.yaml 填写真实邮箱")

    crossref_cfg = cfg.get("crossref") or {}
    if crossref_cfg.get("enabled", True):
        mailto = (crossref_cfg.get("mailto") or "").strip()
        if not mailto or "example.com" in mailto.lower():
            issues.append("CrossRef mailto 未配置或为占位符，请在 config.yaml 填写真实邮箱")

    ai_cfg = cfg.get("ai_screening") or {}
    if ai_cfg.get("enabled", False):
        key_env = ai_cfg.get("api_key_env", "NNAA_FETCH_OPENAI_API_KEY")
        has_key = bool(os.environ.get(key_env, "").strip() or (ai_cfg.get("api_key") or "").strip())
        if not has_key and ai_cfg.get("fail_if_missing_key", False):
            issues.append(f"AI 筛选已启用但未设置 API Key（{key_env} 或 config）")

    for msg in issues:
        logger.warning(f"配置检查: {msg}")
    return issues


def reset_state(state_dir: str):
    path = Path(state_dir) / STATE_FILE
    if path.exists():
        path.unlink()
        print(f"已清除运行状态: {path}")


def _first_day_of_month(d: date) -> date:
    return d.replace(day=1)


def _last_day_of_month(d: date) -> date:
    nxt = _first_day_of_month(d)
    if nxt.month == 12:
        nxt = date(nxt.year + 1, 1, 1)
    else:
        nxt = date(nxt.year, nxt.month + 1, 1)
    return nxt - timedelta(days=1)


def _parse_config_date(value) -> Optional[date]:
    if value is None or str(value).strip() == "":
        return None
    return date.fromisoformat(str(value).strip())


def determine_date_range(
    cfg: dict,
    state: dict,
    source: str,
    cli_date_from: Optional[date],
    cli_date_to: Optional[date],
    search_cfg: dict,
    period_override: Optional[str],
    track: Optional[str] = None,
    batch: Optional[dict] = None,
) -> tuple[date, date]:
    """解析抓取日期窗。默认 explicit：使用 config/CLI 指定的起止日期，全源全轨道一致。"""
    del track, batch  # 日期统一由 search 配置或 CLI 决定，不再按轨道/批次拆分

    mode = (period_override or search_cfg.get("date_range_mode") or "explicit").strip().lower()
    cfg_from = _parse_config_date(search_cfg.get("date_from"))
    cfg_to = _parse_config_date(search_cfg.get("date_to"))
    date_from = cli_date_from or cfg_from
    date_to = cli_date_to or cfg_to or date.today()

    if mode == "explicit":
        if not date_from:
            raise ValueError(
                "search.date_range_mode=explicit 时必须指定 date_from "
                "（config.yaml 的 search.date_from 或命令行 --date-from）"
            )
        if date_from > date_to:
            raise ValueError(f"date_from ({date_from}) 不能晚于 date_to ({date_to})")
        return date_from, date_to

    if cli_date_from:
        return cli_date_from, date_to

    if mode == "from_min_date":
        date_from = date.fromisoformat(str(search_cfg["min_date_from"]))
        return date_from, date_to

    if mode == "rolling":
        days = max(1, int(search_cfg.get("rolling_days", 730)))
        date_from = date_to - timedelta(days=days)
        min_from = search_cfg.get("min_date_from")
        if min_from:
            floor = date.fromisoformat(str(min_from))
            if date_from > floor:
                date_from = floor
        return date_from, date_to

    if mode == "this_month":
        return _first_day_of_month(date_to), date_to

    if mode in ("prev_month", "previous_month", "last_month"):
        first_this = _first_day_of_month(date_to)
        last_prev = first_this - timedelta(days=1)
        first_prev = _first_day_of_month(last_prev)
        return first_prev, last_prev

    last_run_str = state.get(f"{source}_last_run")
    if last_run_str:
        try:
            date_from = date.fromisoformat(last_run_str)
            if date_from > date_to:
                date_from = date_to
            return date_from, date_to
        except ValueError:
            pass

    lookback = int(cfg.get(source, {}).get("initial_lookback_days", 30))
    return date_to - timedelta(days=lookback), date_to


def should_advance_run_state(search_cfg: dict, period_override: Optional[str]) -> bool:
    mode = (period_override or search_cfg.get("date_range_mode") or "explicit").strip().lower()
    return mode == "incremental"


def effective_date_mode(search_cfg: dict, period_override: Optional[str]) -> str:
    return (period_override or search_cfg.get("date_range_mode") or "explicit").strip().lower()


def _partition_records_by_track(records: List[Dict]) -> Dict[str, List[Dict]]:
    """按 nnaa_track 分桶；未知轨道归入 synthesis 合并桶。"""
    buckets: Dict[str, List[Dict]] = {track: [] for track in ALL_NNAA_TRACKS}
    for rec in records:
        track = str(rec.get("nnaa_track") or rec.get("enzyme_track") or "").strip().lower()
        if track in buckets:
            buckets[track].append(rec)
            continue
        cat = str(rec.get("ai_nnaa_category") or rec.get("ai_enzyme_category") or "").strip().lower()
        if cat in buckets:
            buckets[cat].append(rec)
        else:
            buckets.setdefault("pathway", []).append(rec)
    return buckets


def resolve_run_output_dir(output_dir: str, output_cfg: dict, run_time: datetime) -> Path:
    pattern = output_cfg.get("run_dir_pattern", "%Y-%m-%d")
    base = Path(output_dir)
    base.mkdir(parents=True, exist_ok=True)
    name = run_time.strftime(pattern)
    run_dir = base / name
    if run_dir.exists() and any(run_dir.glob("*.csv")):
        name = run_time.strftime(f"{pattern}_%H%M%S")
        run_dir = base / name
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def _write_csv(records: List[Dict], columns: List[str], path: Path, run_time: datetime) -> None:
    rows = []
    for rec in records:
        row = sanitize_record_for_csv({col: rec.get(col, "") for col in columns})
        row["fetched_at"] = run_time.strftime("%Y-%m-%d %H:%M:%S")
        rows.append(row)
    out_cols = columns if "fetched_at" in columns else columns + ["fetched_at"]
    df = pd.DataFrame(rows, columns=out_cols)
    df.to_csv(
        path,
        index=False,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )


def save_output(records: List[Dict], output_cfg: dict, output_dir: str) -> Optional[str]:
    if not records:
        if not output_cfg.get("save_empty_result", False):
            return None
        records = []

    columns = list(
        output_cfg.get(
            "columns",
            ["doi", "title", "source", "pmid", "journal", "pub_date"],
        )
    )
    run_time = datetime.now()
    run_dir = resolve_run_output_dir(output_dir, output_cfg, run_time)
    track_filenames = output_cfg.get("track_filenames") or {
        "pathway": "pathway.csv",
        "enzymatic": "synthesis_enzymatic.csv",
        "fermentation": "synthesis_fermentation.csv",
        "chemical": "synthesis_chemical.csv",
        "hybrid": "synthesis_hybrid.csv",
    }

    track_buckets = _partition_records_by_track(records)
    if not any(track_buckets.values()) and not output_cfg.get("save_empty_result", False):
        return None

    paths_written = []
    for track, track_records in track_buckets.items():
        fname = track_filenames.get(track) or f"{track}.csv"
        track_path = run_dir / fname
        _write_csv(track_records, columns, track_path, run_time)
        paths_written.append(track_path)

    if output_cfg.get("save_combined_in_run_dir", False):
        combined_path = run_dir / output_cfg.get("combined_filename", "all.csv")
        _write_csv(records, columns, combined_path, run_time)
        paths_written.append(combined_path)

    if not paths_written:
        return None
    if output_cfg.get("save_legacy_flat_file", False):
        pattern = output_cfg.get("filename_pattern", "new_dois_%Y-%m_%H%M%S.csv")
        legacy_path = Path(output_dir) / run_time.strftime(pattern)
        _write_csv(records, columns, legacy_path, run_time)

    return str(run_dir)


def run(args):
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"配置文件不存在: {config_path}", file=sys.stderr)
        sys.exit(1)

    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    paths_cfg  = cfg.get("paths", {})
    search_cfg = cfg.get("search", {})
    dedup_cfg  = cfg.get("deduplication", {})
    output_cfg = cfg.get("output", {})
    ai_cfg     = cfg.get("ai_screening", {})

    log_dir    = paths_cfg.get("log_dir", str(ROOT_DIR / "logs"))
    state_dir  = paths_cfg.get("state_dir", str(ROOT_DIR / "state"))
    output_dir = paths_cfg.get("output_dir", str(ROOT_DIR / "output"))
    master_csv = paths_cfg.get("master_csv", "")
    master_doi_col = paths_cfg.get("master_doi_column", "DOI")
    extra_doi_csvs = paths_cfg.get("extra_doi_csvs", [])

    logger = setup_logging(log_dir, verbose=args.verbose)
    logger.info("=" * 60)
    logger.info("非天然氨基酸文献自动获取任务启动")
    logger.info(f"配置文件: {config_path}")

    if args.reset_state:
        reset_state(state_dir)

    if args.dry_run:
        logger.info("[DRY RUN 模式] 仅显示配置，不实际请求")
        validate_config(cfg, logger)
        logger.info(f"search.date_range_mode: {search_cfg.get('date_range_mode', 'explicit')}")
        logger.info(
            f"search.date_from/date_to: {search_cfg.get('date_from')} ~ {search_cfg.get('date_to')}"
        )
        logger.info(f"命令行 --date-from/--date-to: {getattr(args, 'date_from', None)} ~ {getattr(args, 'date_to', None)}")
        logger.info(f"命令行 --period: {getattr(args, 'period', None)}")
        logger.info(f"命令行 --source: {getattr(args, 'source', 'all')}")
        logger.info(f"命令行 --nnaa-track: {getattr(args, 'nnaa_track', 'all')}")
        state = load_state(state_dir)
        period_ov = getattr(args, "period", None)
        cli_date_from = date.fromisoformat(args.date_from) if args.date_from else None
        cli_date_to = date.fromisoformat(args.date_to) if args.date_to else None
        try:
            run_from, run_to = determine_date_range(
                cfg, state, "crossref", cli_date_from, cli_date_to, search_cfg, period_ov
            )
            logger.info(f"预计抓取日期范围: {run_from} ~ {run_to}（PubMed / CrossRef 全轨道共用）")
        except ValueError as e:
            logger.error(str(e))
            return
        nnaa_tracks = resolve_nnaa_tracks(getattr(args, "nnaa_track", "all"))
        if cfg.get("pubmed", {}).get("enabled", True):
            for track in nnaa_tracks:
                pm_q = resolve_pubmed_queries(cfg, track)
                if not pm_q:
                    continue
                logger.info(f"PubMed[{track}] 查询 ({len(pm_q)}): {pm_q}")
        for track in nnaa_tracks:
            crossref_batches = resolve_crossref_query_batches(cfg, track, logger)
            if not crossref_batches:
                continue
            interval = float((cfg.get("crossref") or {}).get("request_interval_sec") or 1.0)
            total_q = sum(len(b.get("queries") or []) for b in crossref_batches)
            for batch in crossref_batches:
                batch_queries = batch.get("queries") or []
                if not batch_queries:
                    continue
                batch_name = batch.get("name") or "batch"
                est_sec = len(batch_queries) * interval
                logger.info(
                    f"CrossRef[{track}/{batch_name}] 检索式 {len(batch_queries)} 条，"
                    f"query_param={batch.get('query_param') or 'default'}，"
                    f"按 interval={interval}s 估算约 {est_sec / 60:.1f} 分钟"
                )
                logger.info(f"CrossRef[{track}/{batch_name}] 前 3 条示例: {batch_queries[:3]!r}")
            logger.info(f"CrossRef[{track}] 合计检索式 {total_q} 条（{len(crossref_batches)} 个批次）")
        logger.info(f"主数据文件（local_pdf）: {master_csv}")
        logger.info(f"附加 DOI 文件: {extra_doi_csvs}")
        logger.info(
            f"主库去重: use_master_csv={dedup_cfg.get('use_master_csv', True)}"
        )
        logger.info(
            f"历史 output 去重: include_previous_outputs="
            f"{dedup_cfg.get('include_previous_outputs', True)}, "
            f"glob={dedup_cfg.get('previous_output_glob', 'new_dois_*.csv')!r}"
        )
        logger.info(f"AI 筛选启用: {ai_cfg.get('enabled', False)}")
        tf_cfg = cfg.get("title_filter") or {}
        logger.info(
            f"标题负向词过滤: enabled={tf_cfg.get('enabled', False)}, "
            f"only_keep_passed={tf_cfg.get('only_keep_passed', True)}"
        )
        at_cfg = cfg.get("article_type_filter") or {}
        logger.info(
            f"文献类型过滤: enabled={at_cfg.get('enabled', False)}, "
            f"exclude_patterns={len(at_cfg.get('exclude_publication_type_patterns') or []) or 'default'}"
        )
        lc_cfg = cfg.get("literature_cleaner") or {}
        logger.info(
            f"文献清洗(ultimate_final_cleaner): enabled={lc_cfg.get('enabled', False)}, "
            f"track_overrides={list((lc_cfg.get('track_overrides') or {}).keys()) or 'none'}"
        )
        logger.info(f"本地 PDF 标记: enabled={(cfg.get('local_pdf') or {}).get('enabled', True)}")
        logger.info(f"输出目录: {output_dir}")
        return

    issues = validate_config(cfg, logger)
    if any("不存在" in i for i in issues):
        logger.error("关键路径缺失，请修正 config.yaml 后重试")
        sys.exit(1)

    state = load_state(state_dir)
    logger.info(f"上次 PubMed 运行时间: {state.get('pubmed_last_run', '首次运行')}")
    logger.info(f"上次 CrossRef 运行时间: {state.get('crossref_last_run', '首次运行')}")

    ai_only = getattr(args, "ai_only", False)
    resume_ai = getattr(args, "resume_ai", False)
    clear_ai_checkpoint_flag = getattr(args, "clear_ai_checkpoint", False)

    if ai_only:
        from ai_checkpoint import load_pending

        logger.info("AI-only 模式：跳过抓取，从 state/pending_ai.json 加载去重后记录")
        try:
            cleaned_records, pending_meta = load_pending(state_dir)
        except (FileNotFoundError, ValueError) as e:
            logger.error(str(e))
            sys.exit(1)
        stats = dict(pending_meta.get("dedup_stats") or {})
        title_filter_stats = dict(pending_meta.get("title_filter_stats") or {})
        article_type_stats = dict(pending_meta.get("article_type_stats") or {})
        literature_cleaner_stats = dict(pending_meta.get("literature_cleaner_stats") or {})
        state = dict(pending_meta.get("state_snapshot") or state)
        all_new_records = []
        goto_post_fetch = True
    else:
        goto_post_fetch = False
        article_type_stats = {}

    if not goto_post_fetch:
        cli_date_from = date.fromisoformat(args.date_from) if args.date_from else None
        cli_date_to   = date.fromisoformat(args.date_to)   if args.date_to   else None
        period_ov = getattr(args, "period", None)
        try:
            run_from, run_to = determine_date_range(
                cfg, state, "crossref", cli_date_from, cli_date_to, search_cfg, period_ov
            )
            logger.info(f"本次抓取日期范围: {run_from} ~ {run_to}（PubMed / CrossRef 全轨道共用）")
        except ValueError as e:
            logger.error(str(e))
            sys.exit(1)

        existing_dois = load_existing_dois_from_config(paths_cfg, dedup_cfg)

        all_new_records: List[Dict] = []

        source = (getattr(args, "source", None) or "all").strip().lower()
        resume_crossref = getattr(args, "resume_crossref", False)
        date_mode = effective_date_mode(search_cfg, period_ov)
        nnaa_tracks = resolve_nnaa_tracks(getattr(args, "nnaa_track", "all"))

        run_pubmed = source in ("all", "pubmed")
        run_crossref = source in ("all", "crossref")
        article_type_stats: Dict[str, int] = {}
        literature_cleaner_stats: Dict[str, int] = {}

        for track in nnaa_tracks:
            track_pubmed_queries = resolve_pubmed_queries(cfg, track)
            if (
                track_pubmed_queries
                and cfg.get("pubmed", {}).get("enabled", True)
                and run_pubmed
            ):
                date_from, date_to = determine_date_range(
                    cfg, state, "pubmed", cli_date_from, cli_date_to, search_cfg, period_ov
                )
                logger.info(f"PubMed[{track}] 搜索日期范围: {date_from} ~ {date_to}")

                fetcher = PubMedFetcher(cfg)
                pubmed_records = fetcher.fetch_since(track_pubmed_queries, date_from, date_to)
                for rec in pubmed_records:
                    rec["nnaa_track"] = track
                logger.info(f"PubMed[{track}] 获取原始记录: {len(pubmed_records)}")
                all_new_records.extend(pubmed_records)
                if date_mode == "incremental" and should_advance_run_state(search_cfg, period_ov):
                    state["pubmed_last_run"] = (date_to + timedelta(days=1)).isoformat()
            elif run_pubmed:
                logger.info(f"PubMed[{track}] 无查询词，跳过")

            crossref_batches = resolve_crossref_query_batches(cfg, track, logger)
            if (
                crossref_batches
                and cfg.get("crossref", {}).get("enabled", True)
                and run_crossref
            ):
                crossref_records: List[Dict] = []
                for batch_idx, batch in enumerate(crossref_batches, 1):
                    batch_queries = batch.get("queries") or []
                    if not batch_queries:
                        continue
                    batch_name = batch.get("name") or f"batch{batch_idx}"
                    date_from, date_to = determine_date_range(
                        cfg,
                        state,
                        "crossref",
                        cli_date_from,
                        cli_date_to,
                        search_cfg,
                        period_ov,
                        track=track,
                        batch=batch,
                    )
                    logger.info(
                        f"CrossRef[{track}/{batch_name}] 搜索日期范围: {date_from} ~ {date_to}, "
                        f"检索式 {len(batch_queries)} 条, query_param={batch.get('query_param') or 'default'}"
                    )

                    if getattr(args, "clear_crossref_checkpoint", False) and track == nnaa_tracks[0] and batch_idx == 1:
                        from crossref_checkpoint import clear_checkpoint

                        clear_checkpoint(state_dir)

                    fetcher = CrossRefFetcher({"crossref": merged_crossref_cfg(cfg, track, batch)})
                    batch_records = fetcher.fetch_since(
                        batch_queries,
                        date_from,
                        date_to,
                        state_dir=state_dir,
                        resume=resume_crossref and track == nnaa_tracks[-1] and batch_idx == len(crossref_batches),
                    )
                    crossref_records.extend(batch_records)
                for rec in crossref_records:
                    rec["nnaa_track"] = track
                logger.info(f"CrossRef[{track}] 获取原始记录: {len(crossref_records)}")
                all_new_records.extend(crossref_records)
                if date_mode == "incremental" and should_advance_run_state(search_cfg, period_ov):
                    state["crossref_last_run"] = (date_to + timedelta(days=1)).isoformat()
            elif run_crossref:
                logger.info(f"CrossRef[{track}] 无查询词，跳过")

        if not run_pubmed:
            logger.info("PubMed 已禁用或 --source 跳过，跳过")
        if not run_crossref:
            logger.info("CrossRef 已禁用或 --source 跳过，跳过")

        logger.info(f"合并后原始记录总数: {len(all_new_records)}")
        from article_type_filter import filter_records as filter_article_types
        from literature_cleaner import filter_records as filter_literature_cleaner

        all_new_records, article_type_stats = filter_article_types(all_new_records, cfg)
        all_new_records, literature_cleaner_stats = filter_literature_cleaner(all_new_records, cfg)
        cleaned_records, stats = deduplicate_new_records(all_new_records, existing_dois, dedup_cfg)

        title_filter_stats = {}
        if cleaned_records:
            from title_filter import filter_records

            cleaned_records, title_filter_stats = filter_records(cleaned_records, cfg)

        local_pdf_dois = set()
        if master_csv and (cfg.get("local_pdf") or {}).get("enabled", True):
            from local_pdf import load_local_pdf_dois, enrich_records as enrich_local_pdf

            local_pdf_dois = load_local_pdf_dois(master_csv, master_doi_col)
            if cleaned_records:
                enrich_local_pdf(cleaned_records, local_pdf_dois)
        elif cleaned_records:
            for rec in cleaned_records:
                rec.setdefault("local_has_pdf", "")

        # 化合物名标注：扫描 title/abstract，填写 nnaa_compound 字段
        if cleaned_records:
            from compound_tagger import enrich_records as enrich_compounds
            enrich_compounds(cleaned_records)
            tagged = sum(1 for r in cleaned_records if r.get("nnaa_compound"))
            logger.info(f"化合物标注完成: {tagged}/{len(cleaned_records)} 条命中至少一个化合物名")

        if cleaned_records:
            from ai_checkpoint import save_pending

            save_pending(
                state_dir,
                cleaned_records,
                meta={
                    "dedup_stats": stats,
                    "title_filter_stats": title_filter_stats,
                    "article_type_stats": article_type_stats,
                    "literature_cleaner_stats": literature_cleaner_stats,
                    "nnaa_tracks": nnaa_tracks,
                    "state_snapshot": state,
                    "date_from": args.date_from or "",
                    "date_to": args.date_to or "",
                },
            )

    ai_stats = {}
    if cleaned_records and ai_cfg.get("enabled", False):
        logger.info("开始 AI 文献筛选…")
        cleaned_records, ai_stats = screen_records(
            cleaned_records,
            ai_cfg,
            state_dir=state_dir,
            resume=resume_ai,
            clear_checkpoint=clear_ai_checkpoint_flag,
        )
        logger.info(f"AI 筛选完成: {ai_stats}")
    elif cleaned_records and not ai_cfg.get("enabled", False):
        for k in ("ai_article_type", "ai_domain_relevant", "ai_has_experiment", "ai_pass", "ai_rationale_zh", "ai_nnaa_category", "ai_synthesis_method", "ai_has_pathway_or_synthesis"):
            for rec in cleaned_records:
                rec.setdefault(k, "")

    if cleaned_records:
        eff_output = dict(output_cfg)
        out_cols = list(output_cfg.get("columns", []))
        for c in (
            "doi",
            "title",
            "source",
            "nnaa_track",
            "pmid",
            "journal",
            "pub_date",
            "publication_types",
            "abstract",
        ):
            if c not in out_cols:
                out_cols.append(c)
        if ai_cfg.get("enabled", False):
            for c in (
                "ai_article_type",
                "ai_domain_relevant",
                "ai_has_experiment",
                "ai_pass",
                "ai_rationale_zh",
                "lit_type_zh",
                "ai_nnaa_category",
                "ai_synthesis_method",
                "ai_has_pathway_or_synthesis",
                "ai_compound",
            ):
                if c not in out_cols:
                    out_cols.append(c)
        # nnaa_compound 始终输出（文本匹配结果，不依赖 AI）
        if "nnaa_compound" not in out_cols:
            out_cols.append("nnaa_compound")
        tf_cfg = cfg.get("title_filter") or {}
        if tf_cfg.get("enabled", False):
            for c in ("title_filter_hits", "title_filter_pass"):
                if c not in out_cols:
                    out_cols.append(c)
        if (cfg.get("local_pdf") or {}).get("enabled", True):
            if "local_has_pdf" not in out_cols:
                out_cols.append("local_has_pdf")
        if "fetched_at" not in out_cols:
            out_cols.append("fetched_at")
        eff_output["columns"] = out_cols
        output_path = save_output(cleaned_records, eff_output, output_dir)
        if output_path:
            track_counts = {
                t: sum(1 for r in cleaned_records if str(r.get("nnaa_track") or "").lower() == t)
                for t in ALL_NNAA_TRACKS
            }
            logger.info(
                f"新增文献已保存至目录: {output_path}  (共 {len(cleaned_records)} 条; "
                f"轨道分布={track_counts})"
            )
        else:
            logger.info("无新增文献，未生成输出文件")
    else:
        logger.info("本次无新增文献（已全部去重）")
        if output_cfg.get("save_empty_result", False):
            save_output([], output_cfg, output_dir)

    save_state(state_dir, state)

    logger.info("任务完成")
    logger.info(
        f"统计: {stats}"
        + (f", 文献类型过滤: {article_type_stats}" if article_type_stats.get("article_type_dropped") else "")
        + (f", 文献清洗: {literature_cleaner_stats}" if literature_cleaner_stats.get("literature_cleaner_dropped") else "")
        + (f", 标题过滤: {title_filter_stats}" if title_filter_stats else "")
        + (f", AI: {ai_stats}" if ai_stats else "")
    )
    logger.info("=" * 60)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=str(ROOT_DIR / "config.yaml"))
    parser.add_argument("--reset-state", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--date-from", metavar="YYYY-MM-DD", default=None,
                        help="抓取起始日期（覆盖 config.yaml 的 search.date_from）")
    parser.add_argument("--date-to", metavar="YYYY-MM-DD", default=None,
                        help="抓取截止日期（覆盖 config.yaml 的 search.date_to，默认今天）")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(
        "--period",
        choices=["incremental", "rolling", "this_month", "prev_month"],
        default=None,
    )
    parser.add_argument(
        "--nnaa-track",
        choices=["pathway", "enzymatic", "fermentation", "chemical", "hybrid", "all"],
        default="all",
        help="检索轨道：pathway=代谢通路, enzymatic=酶法, fermentation=发酵, chemical=化学, hybrid=联用, all=全部（默认 all）",
    )
    parser.add_argument(
        "--source",
        choices=["all", "pubmed", "crossref"],
        default="all",
        help="仅运行指定数据源（默认 all）",
    )
    parser.add_argument(
        "--resume-crossref",
        action="store_true",
        help="从 CrossRef checkpoint 续跑未完成的检索式",
    )
    parser.add_argument(
        "--clear-crossref-checkpoint",
        action="store_true",
        help="清除 CrossRef checkpoint 后重新跑（与 --resume-crossref 互斥）",
    )
    parser.add_argument(
        "--ai-only",
        action="store_true",
        help="跳过抓取/去重，从 state/pending_ai.json 直接跑 AI 及后续 enrich",
    )
    parser.add_argument(
        "--resume-ai",
        action="store_true",
        help="从 state/ai_checkpoint.json 续跑未完成的 AI 批次",
    )
    parser.add_argument(
        "--clear-ai-checkpoint",
        action="store_true",
        help="清除 AI checkpoint 后重新筛（pending_ai.json 仍保留）",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())

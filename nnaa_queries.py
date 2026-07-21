from __future__ import annotations

import logging
from typing import List, Optional

logger = logging.getLogger(__name__)

VALID_NNAA_TRACKS = frozenset({"pathway", "enzymatic", "fermentation", "chemical", "hybrid", "gce", "all"})
ALL_NNAA_TRACKS = ["pathway", "enzymatic", "fermentation", "chemical", "hybrid", "gce"]


def resolve_nnaa_tracks(cli_value: Optional[str]) -> List[str]:
    value = (cli_value or "all").strip().lower()
    if value == "all":
        return list(ALL_NNAA_TRACKS)
    if value in ALL_NNAA_TRACKS:
        return [value]
    raise ValueError(
        f"无效的 --nnaa-track: {cli_value!r}，可选 "
        "pathway | enzymatic | fermentation | chemical | hybrid | gce | all"
    )


def _pubmed_cfg(cfg: dict) -> dict:
    return cfg.get("pubmed") or {}


def _has_split_pubmed(cfg: dict) -> bool:
    pm = _pubmed_cfg(cfg)
    return any(pm.get(f"{track}_queries") for track in ALL_NNAA_TRACKS)


def resolve_pubmed_queries(cfg: dict, nnaa_track: str) -> List[str]:
    pm = _pubmed_cfg(cfg)
    use_code = pm.get("use_keywords_from_code", True)

    key = f"{nnaa_track}_queries"
    queries = [q for q in (pm.get(key) or []) if str(q).strip()]
    if queries:
        return queries

    if use_code:
        from nnaa_keywords import get_pubmed_queries

        code_q = get_pubmed_queries(nnaa_track)
        if code_q:
            return code_q

    legacy = [q for q in (pm.get("queries") or []) if str(q).strip()]
    if legacy and not _has_split_pubmed(cfg):
        return legacy
    return []


def _crossref_base_cfg(cfg: dict) -> dict:
    cr = cfg.get("crossref") or {}
    return {k: v for k, v in cr.items() if k not in ALL_NNAA_TRACKS}


def _crossref_track_cfg(cfg: dict, nnaa_track: str) -> dict:
    cr = cfg.get("crossref") or {}
    track_cfg = cr.get(nnaa_track)
    if isinstance(track_cfg, dict) and track_cfg:
        merged = _crossref_base_cfg(cfg)
        merged.update(track_cfg)
        return merged

    merged = dict(_crossref_base_cfg(cfg))
    merged["query_source"] = merged.get("query_source") or "nnaa_strategy"
    return merged


def _optional_int(value) -> Optional[int]:
    if value is None or str(value).strip() == "":
        return None
    return int(value)


def resolve_crossref_queries(cfg: dict, nnaa_track: str, log: Optional[logging.Logger] = None) -> List[str]:
    log = log or logger
    cr = _crossref_track_cfg(cfg, nnaa_track)
    src = (cr.get("query_source") or "nnaa_strategy").strip().lower()

    if src in ("nnaa_strategy", "nnaa", "default"):
        from nnaa_query_strategy import TRACK_QUERY_BUILDERS

        sub = cr.get("nnaa_strategy") or {}
        max_queries = _optional_int(sub.get("max_queries"))
        builder = TRACK_QUERY_BUILDERS.get(nnaa_track)
        if not builder:
            return []
        qs = builder(max_queries=max_queries)
        log.info(
            "CrossRef[%s] 使用 nnaa_strategy，共 %d 条检索式 (max_queries=%r)",
            nnaa_track,
            len(qs),
            max_queries,
        )
        return qs

    manual = [q for q in (cr.get("queries") or []) if str(q).strip()]
    log.info("CrossRef[%s] 使用 manual 检索式 %d 条", nnaa_track, len(manual))
    return manual


def resolve_crossref_query_batches(
    cfg: dict,
    nnaa_track: str,
    log: Optional[logging.Logger] = None,
) -> List[dict]:
    log = log or logger
    # 支持轨道级别 enabled: false 关闭 CrossRef
    track_raw = (cfg.get("crossref") or {}).get(nnaa_track)
    if isinstance(track_raw, dict) and track_raw.get("enabled") is False:
        log.info("CrossRef[%s] 轨道已禁用（enabled: false），跳过", nnaa_track)
        return []
    cr = _crossref_track_cfg(cfg, nnaa_track)
    batches: List[dict] = []

    main_queries = resolve_crossref_queries(cfg, nnaa_track, log)
    if main_queries:
        batches.append(
            {
                "name": (cr.get("query_source") or "nnaa_strategy"),
                "queries": main_queries,
                "query_param": cr.get("query_param"),
                "limit_per_query": cr.get("limit_per_query"),
            }
        )
    return batches


def merged_crossref_cfg(cfg: dict, nnaa_track: str, batch: Optional[dict] = None) -> dict:
    base = dict(cfg.get("crossref") or {})
    track = _crossref_track_cfg(cfg, nnaa_track)
    merged = dict(base)
    merged.update(track)
    if batch:
        if batch.get("query_param"):
            merged["query_param"] = batch["query_param"]
        if batch.get("limit_per_query") is not None:
            merged["limit_per_query"] = batch["limit_per_query"]
    return merged

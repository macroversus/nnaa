from __future__ import annotations

import logging
import re
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

DEFAULT_EXCLUDE_PATTERNS = [
    r"\breview\b",
    r"systematic review",
    r"meta-?analysis",
    r"\beditorial\b",
    r"\bcomment\b",
    r"\bletter\b",
    r"published erratum",
    r"retraction of publication",
    r"retracted publication",
]


def _compile_patterns(patterns: List[str]) -> List[re.Pattern]:
    out: List[re.Pattern] = []
    for p in patterns:
        p = (p or "").strip()
        if p:
            out.append(re.compile(p, re.IGNORECASE))
    return out


def publication_types_excluded(publication_types: str, patterns: List[re.Pattern]) -> bool:
    text = (publication_types or "").strip()
    if not text:
        return False
    return any(p.search(text) for p in patterns)


def filter_records(
    records: List[Dict],
    cfg: dict,
) -> Tuple[List[Dict], Dict[str, int]]:
    filt_cfg = cfg.get("article_type_filter") or {}
    stats = {
        "article_type_input": len(records),
        "article_type_dropped": 0,
        "article_type_kept": 0,
    }
    if not filt_cfg.get("enabled", False):
        stats["article_type_kept"] = len(records)
        return records, stats

    raw_patterns = list(filt_cfg.get("exclude_publication_type_patterns") or DEFAULT_EXCLUDE_PATTERNS)
    patterns = _compile_patterns(raw_patterns)
    kept: List[Dict] = []

    for rec in records:
        pt = rec.get("publication_types") or ""
        if publication_types_excluded(pt, patterns):
            stats["article_type_dropped"] += 1
            continue
        kept.append(rec)

    stats["article_type_kept"] = len(kept)
    if stats["article_type_dropped"]:
        logger.info(
            "文献类型过滤: 输入=%d, 剔除=%d, 保留=%d",
            stats["article_type_input"],
            stats["article_type_dropped"],
            stats["article_type_kept"],
        )
    return kept, stats

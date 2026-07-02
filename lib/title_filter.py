from __future__ import annotations

import logging
import re
from typing import List, Tuple

logger = logging.getLogger(__name__)

DEFAULT_NEGATIVE_KEYWORDS = [
    "patient",
    "clinical trial",
    "cohort",
    "prognosis",
    "diagnosis",
    "epidemiology",
    "public health",
    "nursing",
    "hospital",
    "systematic review",
    "meta-analysis",
    "economic",
    "social",
    "policy",
    "market",
    "governance",
    "geological",
    "ecological",
    "soil",
    "oceanography",
    "fossil",
    "environmental monitoring",
    "pollution",
    "alloy",
    "ceramic",
    "composite material",
    "semiconductor",
    "quantum",
    "magnetic field",
    "algorithm",
    "software",
]


DEFAULT_POSITIVE_KEYWORDS = [
    "enzyme",
    "enzymatic",
    "biocatalysis",
    "nanozyme",
    "nanozymes",
    "peroxidase-like",
    "oxidase-like",
    "artificial enzyme",
    "enzyme mimic",
    "microreactor",
    "bioreactor",
    "catalysis",
    "catalytic",
]


def _find_hits(text: str, keywords: List[str]) -> List[str]:
    if not text:
        return []
    padded = f" {text.lower()} "
    hits = []
    for kw in keywords:
        k = kw.lower().strip()
        if not k:
            continue
        if f" {k} " in padded:
            hits.append(kw)
    return hits


def filter_records(records: List[dict], cfg: dict) -> Tuple[List[dict], dict]:
    tf = cfg.get("title_filter") or {}
    stats = {
        "input": len(records),
        "dropped": 0,
        "kept": 0,
        "tagged": 0,
        "positive_miss": 0,
    }

    field = tf.get("match_field") or "title"
    negative_kw = list(tf.get("keywords") or DEFAULT_NEGATIVE_KEYWORDS)
    positive_kw = list(tf.get("positive_keywords") or [])
    positive_enabled = tf.get("require_positive_keywords", False)
    only_keep = tf.get("only_keep_passed", True)
    filter_enabled = tf.get("enabled", False) or positive_enabled

    if not filter_enabled:
        for rec in records:
            rec.setdefault("title_filter_hits", "")
            rec.setdefault("title_filter_pass", "true")
        stats["kept"] = len(records)
        return records, stats

    if positive_enabled and not positive_kw:
        positive_kw = list(DEFAULT_POSITIVE_KEYWORDS)

    kept: List[dict] = []
    for rec in records:
        text = str(rec.get(field) or "")
        neg_hits = _find_hits(text, negative_kw) if tf.get("enabled", False) else []
        pos_hits = _find_hits(text, positive_kw) if positive_enabled else []

        hit_parts = []
        if neg_hits:
            hit_parts.append("neg:" + ", ".join(neg_hits))
        if positive_enabled and not pos_hits:
            hit_parts.append("pos:missing")
        elif pos_hits:
            hit_parts.append("pos:" + ", ".join(pos_hits))

        rec["title_filter_hits"] = "; ".join(hit_parts)
        passed = True
        if neg_hits:
            passed = False
            stats["tagged"] += 1
        if positive_enabled and not pos_hits:
            passed = False
            stats["positive_miss"] += 1

        rec["title_filter_pass"] = "true" if passed else "false"
        if not passed and only_keep:
            stats["dropped"] += 1
            continue
        kept.append(rec)
        stats["kept"] += 1

    logger.info(
        f"标题过滤: 输入={stats['input']}, 负向命中={stats['tagged']}, "
        f"正向未命中={stats['positive_miss']}, 丢弃={stats['dropped']}, 保留={stats['kept']}"
    )
    return kept, stats

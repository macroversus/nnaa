from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

DEFAULT_TITLE_PREFIXES = [
    "(invited)",
    "〈",
    "review for",
    "decision letter for",
    "author response for",
    "correction:",
    "retraction:",
    "additions and corrections",
    "editorial:",
    "author correction",
    "review ",
    "a review ",
    "cover feature:",
    "cover picture:",
    "front cover:",
    "frontispiece:",
    "frontispiz:",
    "back cover:",
    "abstract",
    "cheminform abstract:",
    "china & japan:",
    "china:",
    "us:",
    "china's",
    "brazil:",
    "recent",
    "research progress",
    "withdrawn:",
    "retracted",
    "comment",
    "reply",
]

DEFAULT_CONFERENCE_KEYWORDS = ["conference", "mtgabs"]

DEFAULT_DOI_PATTERNS = [r"\.s00", r"/review", r"/decision", r"/response", r"/pdb"]

DEFAULT_TITLE_REGEX_PATTERNS = [
    r"\breview\b",
    r"\brecent\b",
    r"\boverview\b",
]

DEFAULT_COVER_SUFFIX_PATTERN = r"\(.* \d+/\d{4}\)$"


def _track_cfg(cleaner_cfg: dict, enzyme_track: str) -> dict:
    overrides = cleaner_cfg.get("track_overrides") or {}
    base = {k: v for k, v in cleaner_cfg.items() if k != "track_overrides"}
    track = overrides.get(enzyme_track) or {}
    merged = dict(base)
    merged.update(track)
    return merged


def _pub_year(rec: Dict) -> Optional[int]:
    pub_date = str(rec.get("pub_date") or "").strip()
    if not pub_date:
        return None
    m = re.match(r"(\d{4})", pub_date)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def _drop_reason(rec: Dict, enzyme_track: str, cfg: dict) -> Optional[str]:
    track_cfg = _track_cfg(cfg, enzyme_track)
    title = str(rec.get("title") or "")
    lower_title = title.lower()
    doi = str(rec.get("doi") or "").lower()

    if track_cfg.get("require_journal_article", True):
        pub_type = str(rec.get("publication_types") or "").lower().replace("-", " ")
        if pub_type and "journal article" not in pub_type:
            return "non_journal_article"

    prefixes = list(track_cfg.get("title_prefixes") or DEFAULT_TITLE_PREFIXES)
    if any(lower_title.startswith(p.lower()) for p in prefixes if p):
        return "title_prefix"

    for pat in track_cfg.get("title_regex_patterns") or DEFAULT_TITLE_REGEX_PATTERNS:
        if re.search(pat, lower_title, re.IGNORECASE):
            return "title_regex"

    if track_cfg.get("exclude_nanozyme_in_title", enzyme_track == "bio"):
        if "nanozyme" in lower_title:
            return "nanozyme_in_title"

    conf_kw = track_cfg.get("conference_keywords") or DEFAULT_CONFERENCE_KEYWORDS
    conf_pat = "|".join(re.escape(k) for k in conf_kw if k)
    if conf_pat:
        if re.search(conf_pat, lower_title, re.IGNORECASE) or re.search(conf_pat, doi, re.IGNORECASE):
            return "conference_keyword"

    doi_patterns = track_cfg.get("doi_patterns") or DEFAULT_DOI_PATTERNS
    for pat in doi_patterns:
        if re.search(pat, doi, re.IGNORECASE):
            return "doi_pattern"

    if track_cfg.get("drop_numeric_title_start", True) and re.match(r"^\d", title):
        return "numeric_title_start"

    min_words = int(track_cfg.get("min_title_words", 4))
    if min_words > 0 and len(title.split()) < min_words:
        return "short_title"

    cover_pat = track_cfg.get("cover_suffix_pattern") or DEFAULT_COVER_SUFFIX_PATTERN
    if cover_pat and re.search(cover_pat, lower_title):
        return "cover_suffix"

    min_year = track_cfg.get("min_pub_year")
    if min_year is not None:
        year = _pub_year(rec)
        if year is not None and year < int(min_year):
            return "pub_year_before_min"

    return None


def filter_records(
    records: List[Dict],
    cfg: dict,
) -> Tuple[List[Dict], Dict[str, int]]:
    cleaner_cfg = cfg.get("literature_cleaner") or {}
    stats = {
        "literature_cleaner_input": len(records),
        "literature_cleaner_dropped": 0,
        "literature_cleaner_kept": 0,
    }
    reason_counts: Dict[str, int] = {}

    if not cleaner_cfg.get("enabled", False):
        stats["literature_cleaner_kept"] = len(records)
        return records, stats

    kept: List[Dict] = []
    for rec in records:
        track = str(rec.get("nnaa_track") or rec.get("enzyme_track") or "pathway").strip().lower()
        reason = _drop_reason(rec, track, cleaner_cfg)
        if reason:
            stats["literature_cleaner_dropped"] += 1
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
            continue
        kept.append(rec)

    stats["literature_cleaner_kept"] = len(kept)
    for reason, count in sorted(reason_counts.items()):
        stats[f"literature_cleaner_drop_{reason}"] = count

    if stats["literature_cleaner_dropped"]:
        logger.info(
            "文献清洗(ultimate_final_cleaner): 输入=%d, 剔除=%d, 保留=%d, 原因=%s",
            stats["literature_cleaner_input"],
            stats["literature_cleaner_dropped"],
            stats["literature_cleaner_kept"],
            reason_counts,
        )
    return kept, stats

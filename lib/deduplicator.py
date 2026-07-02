import re
import logging
from pathlib import Path

import pandas as pd
from typing import List, Dict, Set, Tuple, Optional

logger = logging.getLogger(__name__)

_DOI_PATTERN = re.compile(r"^10\.\d{4,9}/\S+$")
_NON_WORD = re.compile(r"[^\w\s]+", re.UNICODE)
_WS = re.compile(r"\s+")


def normalize_title(title: str) -> str:
    if not title or not str(title).strip():
        return ""
    t = str(title).lower().strip()
    t = _NON_WORD.sub(" ", t)
    t = _WS.sub(" ", t).strip()
    return t


def normalize_pub_date(pub_date: str) -> str:
    """归一化发表时间，优先 YYYY-MM，否则保留年份或原文。"""
    if not pub_date or not str(pub_date).strip():
        return ""
    raw = str(pub_date).strip()
    m = re.search(r"(20\d{2}|19\d{2})", raw)
    if not m:
        return raw.lower()
    year = m.group(1)
    month_map = {
        "jan": "01", "feb": "02", "mar": "03", "apr": "04",
        "may": "05", "jun": "06", "jul": "07", "aug": "08",
        "sep": "09", "oct": "10", "nov": "11", "dec": "12",
    }
    lower = raw.lower()
    for abbr, mm in month_map.items():
        if abbr in lower:
            return f"{year}-{mm}"
    m2 = re.search(rf"{year}[-/.](\d{{1,2}})", raw)
    if m2:
        return f"{year}-{int(m2.group(1)):02d}"
    return year


def sanitize_csv_text(val) -> str:
    """Normalize text for CSV export; lone \\r breaks row boundaries without quoting."""
    if val is None:
        return ""
    if isinstance(val, float) and val != val:
        return ""
    s = str(val)
    if s.lower() == "nan":
        return ""
    return s.replace("\r\n", "\n").replace("\r", " ")


def sanitize_record_for_csv(rec: Dict) -> Dict:
    return {k: sanitize_csv_text(v) for k, v in rec.items()}


def normalize_doi(doi: str) -> str:
    if not doi:
        return ""
    doi = doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/",
                   "https://dx.doi.org/", "http://dx.doi.org/",
                   "doi:", "DOI:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]
            break
    return doi.lower().strip()


def is_valid_doi(doi: str) -> bool:
    if not doi or len(doi) < 8:
        return False
    return bool(_DOI_PATTERN.match(doi))


def _resolve_doi_column(csv_path: str, preferred: str) -> Optional[str]:
    try:
        cols = list(pd.read_csv(csv_path, nrows=0).columns)
    except Exception:
        return None
    if preferred in cols:
        return preferred
    for fallback in ("doi", "DOI", "Doi"):
        if fallback in cols:
            return fallback
    return None


def _read_doi_column(csv_path: str, doi_column: str) -> List[str]:
    col = _resolve_doi_column(csv_path, doi_column)
    if not col:
        raise KeyError(f"未找到 DOI 列 {doi_column!r}（可用列: {list(pd.read_csv(csv_path, nrows=0).columns)}）")
    df = pd.read_csv(csv_path, usecols=[col], dtype=str)
    return df[col].dropna().tolist()


def load_existing_dois(master_csv: str, doi_column: str = "DOI") -> Set[str]:
    if not master_csv or not str(master_csv).strip():
        return set()
    logger.info(f"加载已有 DOI 列表: {master_csv}")
    try:
        raw_dois = _read_doi_column(master_csv, doi_column)
    except Exception as e:
        logger.error(f"读取主数据文件失败: {e}")
        return set()

    normalized = {normalize_doi(d) for d in raw_dois if d.strip()}
    normalized.discard("")
    logger.info(f"已有 DOI 总数（规范化后）: {len(normalized)}")
    return normalized


def load_existing_dois_union(
    master_csv: Optional[str],
    master_doi_column: str,
    extra_sources: Optional[List[dict]],
) -> Set[str]:
    merged: Set[str] = set()
    if master_csv:
        merged |= load_existing_dois(master_csv, master_doi_column)
    for src in extra_sources or []:
        p = (src or {}).get("path") or ""
        col = (src or {}).get("doi_column") or "doi"
        if not str(p).strip():
            continue
        logger.info(f"加载附加 DOI 列表: {p} (列: {col})")
        try:
            raw = _read_doi_column(p, col)
        except Exception as e:
            logger.warning(f"读取附加文件失败，跳过: {p}: {e}")
            continue
        part = {normalize_doi(d) for d in raw if d and str(d).strip()}
        part.discard("")
        before = len(merged)
        merged |= part
        logger.info(f"  附加 {len(part)} 条规范化 DOI，合并后累计 {len(merged)}（新增 {len(merged) - before}）")
    return merged


def _collect_output_csv_paths(out_path: Path, glob_pattern: str) -> List[Path]:
    files = sorted(out_path.glob(glob_pattern))
    if not files:
        files = sorted(out_path.glob(f"**/{glob_pattern}"))
    for name in (
        "pathway.csv",
        "synthesis_enzymatic.csv",
        "synthesis_fermentation.csv",
        "synthesis_chemical.csv",
        "synthesis_hybrid.csv",
        "all.csv",
    ):
        files.extend(sorted(out_path.glob(f"**/{name}")))
    # 去重并保持稳定顺序
    seen: set = set()
    unique: List[Path] = []
    for p in files:
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            unique.append(p)
    return unique


def load_dois_from_output_dir(
    output_dir: str,
    glob_pattern: str = "new_dois_*.csv",
    doi_column: str = "doi",
) -> Set[str]:
    """Load normalized DOIs from previous run CSVs under output_dir."""
    merged: Set[str] = set()
    if not output_dir or not str(output_dir).strip():
        return merged

    out_path = Path(output_dir)
    if not out_path.is_dir():
        logger.info(f"输出目录不存在，跳过历史 output 去重: {output_dir}")
        return merged

    files = _collect_output_csv_paths(out_path, glob_pattern)
    if not files:
        logger.info(f"未找到历史 output 文件 ({glob_pattern} 及 bio/nano 子目录)，跳过")
        return merged

    logger.info(f"扫描历史 output 去重: {len(files)} 个文件 ({glob_pattern} + 日期子目录)")
    for csv_path in files:
        logger.info(f"  读取历史 output: {csv_path}")
        try:
            raw = _read_doi_column(str(csv_path), doi_column)
        except Exception as e:
            logger.warning(f"  读取失败，跳过: {csv_path}: {e}")
            continue
        part = {normalize_doi(d) for d in raw if d and str(d).strip()}
        part.discard("")
        before = len(merged)
        merged |= part
        logger.info(
            f"  {len(part)} 条 DOI，本批新增 {len(merged) - before}，"
            f"历史 output 累计 {len(merged)}"
        )
    return merged


def load_existing_dois_from_config(paths_cfg: dict, dedup_cfg: Optional[dict] = None) -> Set[str]:
    """Union of master, extra CSVs, and optionally previous output CSVs."""
    paths_cfg = paths_cfg or {}
    dedup_cfg = dedup_cfg or {}

    merged: Set[str] = set()
    if dedup_cfg.get("use_master_csv", True):
        merged = load_existing_dois_union(
            paths_cfg.get("master_csv", ""),
            paths_cfg.get("master_doi_column", "DOI"),
            paths_cfg.get("extra_doi_csvs") or [],
        )
    elif paths_cfg.get("extra_doi_csvs"):
        merged = load_existing_dois_union("", "DOI", paths_cfg.get("extra_doi_csvs") or [])
    else:
        logger.info("主库去重已关闭（deduplication.use_master_csv=false），跳过酶库 DOI 去重")

    if not dedup_cfg.get("include_previous_outputs", True):
        return merged

    output_dir = paths_cfg.get("output_dir", "")
    glob_pattern = dedup_cfg.get("previous_output_glob", "new_dois_*.csv")
    doi_column = dedup_cfg.get("previous_output_doi_column", "doi")
    part = load_dois_from_output_dir(output_dir, glob_pattern, doi_column)
    if not part:
        return merged

    before = len(merged)
    merged |= part
    logger.info(
        f"并入历史 output 后已知 DOI 累计 {len(merged)} "
        f"（output 新增 {len(merged) - before}）"
    )
    return merged


def deduplicate_new_records(
    records: List[Dict],
    existing_dois: Set[str],
    dedup_config: dict,
) -> Tuple[List[Dict], dict]:
    do_normalize  = dedup_config.get("normalize_doi", True)
    do_strip      = dedup_config.get("strip_doi_prefix", True)
    skip_invalid  = dedup_config.get("skip_invalid_doi", True)

    stats = {
        "total_raw": len(records),
        "skipped_no_doi": 0,
        "skipped_invalid_doi": 0,
        "skipped_duplicate_in_batch": 0,
        "skipped_duplicate_in_master": 0,
        "skipped_duplicate_title": 0,
        "skipped_duplicate_title_date": 0,
        "new_count": 0,
    }

    dedupe_title = dedup_config.get("dedupe_by_title", False)
    dedupe_title_date = dedup_config.get("dedupe_by_title_and_pub_date", False)

    seen_in_batch: Set[str] = set()
    seen_titles: Set[str] = set()
    seen_title_dates: Set[Tuple[str, str]] = set()
    cleaned: List[Dict] = []

    for rec in records:
        doi = rec.get("doi", "")

        if not doi:
            stats["skipped_no_doi"] += 1
            continue

        if do_normalize or do_strip:
            doi = normalize_doi(doi)

        if skip_invalid and not is_valid_doi(doi):
            stats["skipped_invalid_doi"] += 1
            logger.debug(f"无效 DOI 跳过: {doi!r}")
            continue

        if doi in seen_in_batch:
            stats["skipped_duplicate_in_batch"] += 1
            continue

        if doi in existing_dois:
            stats["skipped_duplicate_in_master"] += 1
            continue

        norm_title = normalize_title(rec.get("title") or "")
        norm_date = normalize_pub_date(rec.get("pub_date") or "")

        if dedupe_title_date and norm_title and norm_date:
            td_key = (norm_title, norm_date)
            if td_key in seen_title_dates:
                stats["skipped_duplicate_title_date"] += 1
                continue

        if dedupe_title and norm_title:
            if norm_title in seen_titles:
                stats["skipped_duplicate_title"] += 1
                continue

        seen_in_batch.add(doi)
        if norm_title:
            seen_titles.add(norm_title)
        if dedupe_title_date and norm_title and norm_date:
            seen_title_dates.add((norm_title, norm_date))

        rec = dict(rec)
        rec["doi"] = doi
        cleaned.append(rec)

    stats["new_count"] = len(cleaned)

    logger.info(
        f"去重清洗结果: 原始={stats['total_raw']}, "
        f"无DOI={stats['skipped_no_doi']}, "
        f"格式无效={stats['skipped_invalid_doi']}, "
        f"批次内DOI重复={stats['skipped_duplicate_in_batch']}, "
        f"已在主库={stats['skipped_duplicate_in_master']}, "
        f"标题重复={stats['skipped_duplicate_title']}, "
        f"标题+日期重复={stats['skipped_duplicate_title_date']}, "
        f"新增={stats['new_count']}"
    )
    return cleaned, stats

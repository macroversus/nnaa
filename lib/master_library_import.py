"""
从酶相关主库（allpdfs）中检索非天然氨基酸相关文献，并入 NNAA 输出。

主库 CSV 仅含 DOI/Filename；标题等元数据来自 paths.master_metadata_csv。
匹配依据：标题中的 NNAA 核心词 + 轨道信号词（与 nnaa_keywords 一致）。
"""

from __future__ import annotations

import csv
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import yaml

from deduplicator import normalize_doi, sanitize_record_for_csv
from nnaa_keywords import (
    CHEMICAL_SIGNALS_EN,
    ENZYMATIC_SIGNALS_EN,
    FERMENTATION_SIGNALS_EN,
    HYBRID_SIGNALS_EN,
    NNAA_CORE_EN,
    PATHWAY_SIGNALS_EN,
)

logger = logging.getLogger(__name__)

TRACK_SIGNALS: Dict[str, Tuple[str, ...]] = {
    "pathway": PATHWAY_SIGNALS_EN,
    "enzymatic": ENZYMATIC_SIGNALS_EN,
    "fermentation": FERMENTATION_SIGNALS_EN,
    "chemical": CHEMICAL_SIGNALS_EN,
    "hybrid": HYBRID_SIGNALS_EN,
}

TRACK_PRIORITY = ("hybrid", "enzymatic", "fermentation", "chemical", "pathway")

# 标题中须出现至少一条 NNAA 核心表述（比单独 "amino acid" 更严格）
NNAA_TITLE_MARKERS: Tuple[str, ...] = tuple(
    dict.fromkeys(
        list(NNAA_CORE_EN)
        + [
            "unnatural amino",
            "non-natural amino",
            "noncanonical amino",
            "non-proteinogenic",
            "non-standard amino",
            "orthogonal aminoacyl",
            "pyrrolysyl-tRNA",
            "amber suppression",
            "ncAA",
        ]
    )
)


def _norm(text: str) -> str:
    return f" {(text or '').lower()} "


def title_has_nnaa_marker(title: str) -> bool:
    t = _norm(title)
    return any(m.lower() in t for m in NNAA_TITLE_MARKERS)


def assign_track_from_title(title: str) -> Optional[str]:
    if not title_has_nnaa_marker(title):
        return None
    t = (title or "").lower()
    for track in TRACK_PRIORITY:
        if any(sig.lower() in t for sig in TRACK_SIGNALS[track]):
            return track
    return "pathway"


def _parse_year(val) -> Optional[int]:
    if val is None or (isinstance(val, float) and val != val):
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        return int(float(s[:4]))
    except ValueError:
        return None


def _year_in_range(year: Optional[int], date_from: Optional[str], date_to: Optional[str]) -> bool:
    if year is None:
        return True
    if date_from:
        try:
            y0 = int(date_from[:4])
            if year < y0:
                return False
        except ValueError:
            pass
    if date_to:
        try:
            y1 = int(date_to[:4])
            if year > y1:
                return False
        except ValueError:
            pass
    return True


def load_master_metadata(paths_cfg: dict) -> pd.DataFrame:
    master_csv = (paths_cfg.get("master_csv") or "").strip()
    meta_csv = (paths_cfg.get("master_metadata_csv") or "").strip()
    doi_col = paths_cfg.get("master_doi_column") or "DOI"

    if not master_csv or not Path(master_csv).exists():
        raise FileNotFoundError(f"主库 CSV 不存在: {master_csv}")
    if not meta_csv or not Path(meta_csv).exists():
        raise FileNotFoundError(f"主库元数据 CSV 不存在: {meta_csv}")

    master = pd.read_csv(master_csv, usecols=[doi_col], dtype=str)
    meta = pd.read_csv(meta_csv, dtype=str, low_memory=False)
    if doi_col not in meta.columns and "DOI" in meta.columns:
        doi_col = "DOI"

    merged = master.merge(meta, on=doi_col, how="inner")
    logger.info(
        f"主库元数据: master={len(master)}, metadata={len(meta)}, 合并={len(merged)}"
    )
    return merged


def scan_master_library(
    cfg: dict,
    *,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
    all_years: bool = False,
) -> List[Dict]:
    paths_cfg = cfg.get("paths") or {}
    search_cfg = cfg.get("search") or {}

    if not all_years:
        date_from = date_from or search_cfg.get("date_from")
        date_to = date_to or search_cfg.get("date_to")

    df = load_master_metadata(paths_cfg)
    doi_col = paths_cfg.get("master_doi_column") or "DOI"
    title_col = "title" if "title" in df.columns else None
    if not title_col:
        raise KeyError("主库元数据缺少 title 列")

    records: List[Dict] = []
    skipped_date = 0
    skipped_no_match = 0

    for _, row in df.iterrows():
        title = str(row.get(title_col) or "").strip()
        if not title_has_nnaa_marker(title):
            skipped_no_match += 1
            continue

        year = _parse_year(row.get("year"))
        if not all_years and not _year_in_range(year, date_from, date_to):
            skipped_date += 1
            continue

        track = assign_track_from_title(title)
        if not track:
            continue

        doi_raw = str(row.get(doi_col) or "").strip()
        doi = normalize_doi(doi_raw)
        if not doi:
            continue

        pub_date = str(row.get("year") or "").strip()
        records.append(
            {
                "doi": doi,
                "title": title,
                "source": "master_library",
                "nnaa_track": track,
                "pmid": "",
                "journal": str(row.get("journal") or "").strip(),
                "pub_date": pub_date,
                "publication_types": "",
                "abstract": "",
                "local_has_pdf": "true",
            }
        )

    logger.info(
        f"主库 NNAA 标题匹配: 候选={len(records)}, "
        f"无 NNAA 词={skipped_no_match}, 日期外={skipped_date}"
    )
    return records


def load_dois_from_output(output_dir: str) -> Set[str]:
    out = Path(output_dir)
    if not out.is_dir():
        return set()

    dois: Set[str] = set()
    for pattern in ("nnaa_dois_*.csv", "**/all.csv", "**/pathway.csv", "**/synthesis_*.csv"):
        for csv_path in out.glob(pattern):
            try:
                df = pd.read_csv(csv_path, usecols=["doi"], dtype=str)
            except Exception:
                try:
                    df = pd.read_csv(csv_path, usecols=["DOI"], dtype=str)
                    df = df.rename(columns={"DOI": "doi"})
                except Exception:
                    continue
            for d in df["doi"].dropna():
                nd = normalize_doi(str(d))
                if nd:
                    dois.add(nd)
    return dois


def filter_new_records(records: List[Dict], existing_dois: Set[str]) -> List[Dict]:
    new: List[Dict] = []
    for rec in records:
        doi = normalize_doi(rec.get("doi") or "")
        if doi and doi not in existing_dois:
            new.append(rec)
    logger.info(f"主库导入: 匹配 {len(records)} 条，已有 output 中 {len(records) - len(new)} 条，新增 {len(new)} 条")
    return new


def find_latest_run_csv(output_dir: Path) -> Optional[Path]:
    legacy = sorted(output_dir.glob("nnaa_dois_*.csv"), reverse=True)
    if legacy:
        return legacy[0]
    run_dirs = sorted(
        [d for d in output_dir.iterdir() if d.is_dir() and d.name != "key_results"],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for run_dir in run_dirs:
        combined = run_dir / "all.csv"
        if combined.exists():
            return combined
    return None


def merge_records_into_output(
    base_csv: Path,
    new_records: List[Dict],
    output_dir: Path,
    output_cfg: dict,
) -> Path:
    if not new_records:
        return base_csv

    columns = list(output_cfg.get("columns") or [])
    for c in ("doi", "title", "source", "nnaa_track", "local_has_pdf", "fetched_at"):
        if c not in columns:
            columns.append(c)

    base_df = pd.read_csv(base_csv, dtype=str, keep_default_na=False)
    run_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    rows = []
    for rec in new_records:
        row = sanitize_record_for_csv({col: rec.get(col, "") for col in columns})
        row["fetched_at"] = run_time
        rows.append(row)

    new_df = pd.DataFrame(rows, columns=columns if "fetched_at" in columns else columns + ["fetched_at"])
    combined = pd.concat([base_df, new_df], ignore_index=True)

    pattern = output_cfg.get("filename_pattern", "nnaa_dois_%Y-%m_%H%M%S.csv")
    out_path = output_dir / datetime.now().strftime(pattern)
    combined.to_csv(
        out_path,
        index=False,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )
    logger.info(f"已合并写入: {out_path}（原 {len(base_df)} + 新增 {len(new_records)} = {len(combined)}）")
    return out_path


def run_import_master(
    config_path: Path,
    *,
    all_years: bool = False,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
    skip_ai: bool = False,
    verbose: bool = False,
) -> dict:
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    paths_cfg = cfg.get("paths") or {}
    output_dir = Path(paths_cfg.get("output_dir") or Path(config_path).parent / "output")
    state_dir = paths_cfg.get("state_dir") or str(Path(config_path).parent / "state")
    ai_cfg = cfg.get("ai_screening") or {}
    output_cfg = cfg.get("output") or {}

    candidates = scan_master_library(cfg, date_from=date_from, date_to=date_to, all_years=all_years)
    existing = load_dois_from_output(str(output_dir))
    new_records = filter_new_records(candidates, existing)

    stats = {
        "candidates": len(candidates),
        "already_in_output": len(candidates) - len(new_records),
        "imported": 0,
        "ai_stats": {},
        "output_csv": "",
    }

    if not new_records:
        logger.info("主库无新增 NNAA 文献需导入")
        return stats

    if ai_cfg.get("enabled", False) and not skip_ai:
        from ai_screener import screen_records

        logger.info(f"对 {len(new_records)} 条主库导入记录运行 AI 筛选…")
        new_records, ai_stats = screen_records(new_records, ai_cfg, state_dir=state_dir)
        stats["ai_stats"] = ai_stats
    else:
        for rec in new_records:
            for k in (
                "ai_article_type",
                "ai_domain_relevant",
                "ai_has_experiment",
                "ai_pass",
                "ai_rationale_zh",
                "ai_nnaa_category",
                "ai_synthesis_method",
                "ai_has_pathway_or_synthesis",
            ):
                rec.setdefault(k, "")

    base_csv = find_latest_run_csv(output_dir)
    if not base_csv:
        logger.warning("未找到已有 output CSV，将仅写入主库导入记录")
        pattern = output_cfg.get("filename_pattern", "nnaa_dois_%Y-%m_%H%M%S.csv")
        out_path = output_dir / datetime.now().strftime(pattern)
        columns = list(output_cfg.get("columns") or [])
        for c in ("doi", "title", "source", "nnaa_track", "local_has_pdf", "fetched_at"):
            if c not in columns:
                columns.append(c)
        run_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        rows = []
        for rec in new_records:
            row = sanitize_record_for_csv({col: rec.get(col, "") for col in columns})
            row["fetched_at"] = run_time
            rows.append(row)
        out_cols = columns if "fetched_at" in columns else columns + ["fetched_at"]
        pd.DataFrame(rows, columns=out_cols).to_csv(
            out_path, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
        )
        stats["output_csv"] = str(out_path)
        stats["imported"] = len(new_records)
        return stats

    out_path = merge_records_into_output(base_csv, new_records, output_dir, output_cfg)
    stats["output_csv"] = str(out_path)
    stats["imported"] = len(new_records)
    return stats

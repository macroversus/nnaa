#!/usr/bin/env python3
"""Repair output CSV: restore pending base fields, merge AI/enrich columns, fix \\r breakage."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import bootstrap  # noqa: F401

ROOT = Path(__file__).resolve().parent.parent

from lib.deduplicator import is_valid_doi, normalize_doi, sanitize_record_for_csv

OVERLAY_PREFIXES = ("ai_",)
OVERLAY_EXACT = {
    "local_has_pdf",
    "title_filter_hits",
    "title_filter_pass",
    "fetched_at",
}


def _should_overlay(key: str) -> bool:
    return key.startswith(OVERLAY_PREFIXES) or key in OVERLAY_EXACT


def repair(csv_path: Path, pending_path: Path, output_path: Path | None = None) -> int:
    pending = json.load(open(pending_path, encoding="utf-8"))
    pending_records = pending.get("records") or []
    if not pending_records:
        raise ValueError(f"pending 为空: {pending_path}")

    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        by_doi = {}
        skipped_invalid = 0
        for row in reader:
            doi = normalize_doi(row.get("doi") or "")
            if not doi or not is_valid_doi(doi):
                skipped_invalid += 1
                continue
            by_doi[doi] = row

    out_rows = []
    missing_ai = 0
    for rec in pending_records:
        doi = normalize_doi(rec.get("doi") or "")
        if not doi:
            continue
        merged = dict(rec)
        src = by_doi.get(doi)
        if src:
            for k, v in src.items():
                if _should_overlay(k):
                    merged[k] = v
        else:
            missing_ai += 1
        out_rows.append(sanitize_record_for_csv(merged))

    if missing_ai:
        print(f"警告: {missing_ai} 条 pending 记录在 CSV 中无 AI 覆盖", file=sys.stderr)
    if skipped_invalid:
        print(f"已跳过无效 DOI 行: {skipped_invalid}", file=sys.stderr)

    out_path = output_path or csv_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in out_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    tmp.replace(out_path)
    print(f"已写入 {out_path} ({len(out_rows)} 条)")
    return len(out_rows)


def main():
    parser = argparse.ArgumentParser(description="修复 output CSV（pending 基准 + AI 列合并）")
    parser.add_argument(
        "-i",
        "--input",
        default=str(SCRIPT_DIR / "output/new_dois_2026-05_143600.csv"),
    )
    parser.add_argument(
        "-p",
        "--pending",
        default=str(SCRIPT_DIR / "state/pending_ai.json"),
    )
    parser.add_argument("-o", "--output", help="输出路径（默认覆盖输入）")
    args = parser.parse_args()
    repair(Path(args.input), Path(args.pending), Path(args.output) if args.output else None)


if __name__ == "__main__":
    main()

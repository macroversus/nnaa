#!/usr/bin/env python3
"""
从 AI 筛选后的 CSV 导出「关键结果文件」。

关键结果定义：ai_pass=true 且 ai_domain_relevant=true（或 ai_pass=true 且非 batch_api_error）。
按 nnaa_track / ai_nnaa_category 分轨输出 usable CSV + 汇总 JSON。
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

import bootstrap  # noqa: F401

from nnaa_queries import ALL_NNAA_TRACKS

SCRIPT_DIR = Path(__file__).parent.resolve()

KEY_TRACK_FILES = {
    "pathway": "关键_代谢通路.csv",
    "enzymatic": "关键_酶法合成.csv",
    "fermentation": "关键_生物发酵.csv",
    "chemical": "关键_化学合成.csv",
    "hybrid": "关键_酶化学联用.csv",
}
ALL_KEY_FILE = "关键_全部可用文献.csv"
SUMMARY_FILE = "关键结果汇总.json"


def _is_usable(row: dict) -> bool:
    ai_pass = str(row.get("ai_pass") or "").strip().lower()
    if ai_pass != "true":
        return False
    rationale = str(row.get("ai_rationale_zh") or "").strip()
    if rationale == "batch_api_error_kept":
        return False
    return True


def _resolve_track(row: dict) -> str:
    track = str(row.get("nnaa_track") or "").strip().lower()
    if track in ALL_NNAA_TRACKS:
        return track
    cat = str(row.get("ai_nnaa_category") or "").strip().lower()
    if cat in ALL_NNAA_TRACKS:
        return cat
    return "pathway"


def find_latest_run_csv(output_dir: Path) -> Optional[Path]:
    legacy = sorted(output_dir.glob("nnaa_dois_*.csv"), reverse=True)
    if legacy:
        return legacy[0]
    run_dirs = sorted(
        [d for d in output_dir.iterdir() if d.is_dir()],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for run_dir in run_dirs:
        combined = run_dir / "all.csv"
        if combined.exists():
            return combined
    return None


def export_key_results(
    input_csv: Path,
    output_dir: Path,
    *,
    run_label: Optional[str] = None,
) -> dict:
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    rows = df.to_dict(orient="records")
    usable = [r for r in rows if _is_usable(r)]

    label = run_label or datetime.now().strftime("%Y-%m-%d_%H%M%S")
    key_dir = output_dir / "key_results" / label
    key_dir.mkdir(parents=True, exist_ok=True)

    by_track: Dict[str, List[dict]] = {t: [] for t in ALL_NNAA_TRACKS}
    for row in usable:
        by_track[_resolve_track(row)].append(row)

    columns = list(df.columns)
    stats = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source_csv": str(input_csv),
        "total_input": len(rows),
        "total_usable": len(usable),
        "by_track": {t: len(by_track[t]) for t in ALL_NNAA_TRACKS},
        "output_dir": str(key_dir),
        "files": {},
    }

    for track, fname in KEY_TRACK_FILES.items():
        path = key_dir / fname
        track_rows = by_track[track]
        pd.DataFrame(track_rows, columns=columns if track_rows else columns).to_csv(
            path,
            index=False,
            encoding="utf-8-sig",
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        stats["files"][track] = str(path)

    all_path = key_dir / ALL_KEY_FILE
    pd.DataFrame(usable, columns=columns if usable else columns).to_csv(
        all_path,
        index=False,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )
    stats["files"]["all"] = str(all_path)

    summary_path = key_dir / SUMMARY_FILE
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    stats["files"]["summary"] = str(summary_path)

    return stats


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="导出非天然氨基酸 AI 筛选后的关键结果 CSV")
    parser.add_argument("-i", "--input", help="输入 CSV（默认取 output 下最新 all.csv 或 nnaa_dois_*.csv）")
    parser.add_argument(
        "-o",
        "--output-dir",
        default=str(SCRIPT_DIR / "output"),
        help="输出根目录（关键文件写入 output/key_results/）",
    )
    parser.add_argument("--label", default=None, help="本次关键结果子目录名")
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    if args.input:
        input_csv = Path(args.input)
    else:
        input_csv = find_latest_run_csv(output_dir)
        if not input_csv:
            print("未找到输入 CSV，请用 -i 指定或先运行 nnaa_fetch.py run", file=sys.stderr)
            return 1

    if not input_csv.exists():
        print(f"输入文件不存在: {input_csv}", file=sys.stderr)
        return 1

    stats = export_key_results(input_csv, output_dir, run_label=args.label)
    print(f"关键结果已导出至: {stats['output_dir']}")
    print(f"  输入 {stats['total_input']} 条 → 可用 {stats['total_usable']} 条")
    for track, n in stats["by_track"].items():
        if n:
            print(f"  {track}: {n} 条 → {stats['files'][track]}")
    print(f"  汇总: {stats['files']['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

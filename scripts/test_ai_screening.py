#!/usr/bin/env python3
"""对已有 CSV 小批量试跑 AI 筛选（不改动 state）。"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

import bootstrap  # noqa: F401

ROOT = Path(__file__).resolve().parent.parent

from lib.ai_screener import screen_records


def _has_api_key(ai_cfg: dict) -> bool:
    key_env = ai_cfg.get("api_key_env", "NNAA_FETCH_OPENAI_API_KEY")
    return bool(
        os.environ.get(key_env, "").strip() or (ai_cfg.get("api_key") or "").strip()
    )


def main():
    parser = argparse.ArgumentParser(description="对 CSV 中文献小批量试跑 AI 筛选")
    parser.add_argument("-c", "--config", default=str(ROOT / "config.yaml"))
    parser.add_argument("-i", "--input", required=True, help="输入 CSV")
    parser.add_argument("-o", "--output", help="输出 CSV（默认输入名加 _ai 后缀）")
    parser.add_argument("--limit", type=int, default=6, help="试跑条数（默认 6）")
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="临时关闭 only_keep_passed，保留未通过项便于查看",
    )
    args = parser.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    ai_cfg = dict(cfg.get("ai_screening") or {})
    if not _has_api_key(ai_cfg):
        print(
            "未检测到 API Key。请设置环境变量 "
            f"{ai_cfg.get('api_key_env', 'NNAA_FETCH_OPENAI_API_KEY')} "
            "或在 config.yaml 的 ai_screening.api_key 中填写。",
            file=sys.stderr,
        )
        sys.exit(1)

    ai_cfg["enabled"] = True
    if args.keep_all:
        ai_cfg["only_keep_passed"] = False

    df = pd.read_csv(args.input, dtype=str)
    rows = df.head(max(1, args.limit)).to_dict(orient="records")

    print(f"试跑 {len(rows)} 条 | model={ai_cfg.get('model', 'gpt-4o-mini')} | base={ai_cfg.get('api_base')}")
    screened, stats = screen_records(rows, ai_cfg)
    print("统计:", stats)

    out = args.output
    if not out:
        p = Path(args.input)
        out = str(p.with_name(p.stem + "_ai" + p.suffix))

    pd.DataFrame(screened).to_csv(out, index=False, encoding="utf-8-sig")
    print(f"已写入: {out}（{len(screened)} 条）")

    cols = [
        "doi",
        "title",
        "ai_article_type",
        "ai_domain_relevant",
        "ai_has_experiment",
        "ai_pass",
        "ai_nnaa_category",
        "ai_synthesis_method",
        "ai_has_pathway_or_synthesis",
    ]
    show = [c for c in cols if c in pd.read_csv(out, nrows=0).columns]
    if show:
        preview = pd.read_csv(out, dtype=str)[show]
        print("\n预览:")
        print(preview.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
非天然氨基酸文献流水线 — 独立入口

流程：检索（PubMed + CrossRef）→ AI 判断 → 关键结果输出

用法：
  python3 nnaa_fetch.py run              # 完整三步流水线
  python3 nnaa_fetch.py search           # 仅检索 + 去重
  python3 nnaa_fetch.py ai               # 仅 AI（基于 state/pending_ai.json）
  python3 nnaa_fetch.py import-master   # 从酶库主库检索 NNAA 相关并合并
  python3 nnaa_fetch.py export           # 仅导出关键 CSV
  python3 nnaa_fetch.py keywords         # 打印各轨道检索关键词/检索式
  python3 nnaa_fetch.py run --dry-run    # 预览配置与检索式
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import yaml

import bootstrap  # noqa: F401  — 注入 ROOT_DIR / lib 到 sys.path

SCRIPT_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = SCRIPT_DIR / "config.yaml"


def _load_cfg(config_path: Path) -> dict:
    with open(config_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _write_temp_cfg(cfg: dict, overrides: dict) -> Path:
    merged = copy.deepcopy(cfg)
    for key, val in overrides.items():
        parts = key.split(".")
        node = merged
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        node[parts[-1]] = val
    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".yaml",
        prefix="nnaa_fetch_",
        delete=False,
        encoding="utf-8",
    )
    yaml.dump(merged, tmp, allow_unicode=True, default_flow_style=False)
    tmp.close()
    return Path(tmp.name)


def _run_fetch(extra_argv: List[str], config_path: Path) -> int:
    import bootstrap  # noqa: F401
    from lib.fetch_literature import parse_args, run

    old_argv = sys.argv
    sys.argv = ["fetch_literature.py", "-c", str(config_path)] + extra_argv
    try:
        run(parse_args())
    finally:
        sys.argv = old_argv
    return 0


def cmd_keywords(_args) -> int:
    from nnaa_keywords import PUBMED_QUERIES_BY_TRACK, CROSSREF_SIGNALS_BY_TRACK, NNAA_CORE_EN
    from nnaa_query_strategy import TRACK_QUERY_BUILDERS

    print("=" * 60)
    print("非天然氨基酸检索关键词库（nnaa_keywords.py）")
    print("=" * 60)
    print("\n【核心同义词】")
    for term in NNAA_CORE_EN:
        print(f"  - {term}")

    for track, pm_q in PUBMED_QUERIES_BY_TRACK.items():
        print(f"\n【轨道: {track}】")
        print(f"  CrossRef 信号词 ({len(CROSSREF_SIGNALS_BY_TRACK[track])} 个):")
        for s in CROSSREF_SIGNALS_BY_TRACK[track][:8]:
            print(f"    · {s}")
        if len(CROSSREF_SIGNALS_BY_TRACK[track]) > 8:
            print(f"    · ... 共 {len(CROSSREF_SIGNALS_BY_TRACK[track])} 个")
        builder = TRACK_QUERY_BUILDERS[track]
        cr_n = len(builder(max_queries=48))
        print(f"  CrossRef 检索式: {cr_n} 条（max 48 示例）")
        print(f"  PubMed 检索式 ({len(pm_q)} 条):")
        for i, q in enumerate(pm_q, 1):
            preview = q if len(q) <= 120 else q[:117] + "..."
            print(f"    {i}. {preview}")
    return 0


def cmd_search(args) -> int:
    config_path = Path(args.config)
    extra = ["--source", args.source]
    if args.date_from:
        extra += ["--date-from", args.date_from]
    if args.date_to:
        extra += ["--date-to", args.date_to]
    if args.track and args.track != "all":
        extra += ["--nnaa-track", args.track]
    if args.dry_run:
        extra.append("--dry-run")
    if args.verbose:
        extra.append("-v")
    cfg = _load_cfg(config_path)
    tmp = _write_temp_cfg(cfg, {"ai_screening.enabled": False})
    try:
        return _run_fetch(extra, tmp)
    finally:
        tmp.unlink(missing_ok=True)


def cmd_ai(args) -> int:
    extra = ["--ai-only"]
    if args.resume:
        extra.append("--resume-ai")
    if args.clear_checkpoint:
        extra.append("--clear-ai-checkpoint")
    if args.verbose:
        extra.append("-v")
    return _run_fetch(extra, Path(args.config))


def cmd_import_master(args) -> int:
    import json
    from lib.master_library_import import run_import_master

    config_path = Path(args.config)
    stats = run_import_master(
        config_path,
        all_years=args.all_years,
        date_from=args.date_from,
        date_to=args.date_to,
        skip_ai=args.skip_ai,
        verbose=args.verbose,
    )
    print(json.dumps(stats, indent=2, ensure_ascii=False))

    if stats.get("imported", 0) > 0 and not args.skip_ai:
        print("\n>>> 重新导出关键结果")
        export_args = argparse.Namespace(
            input=stats.get("output_csv") or args.input,
            output_dir=args.output_dir,
            label=args.label,
        )
        return cmd_export(export_args)
    return 0


def cmd_export(args) -> int:
    from export_key_results import export_key_results, find_latest_run_csv

    output_dir = Path(args.output_dir)
    if args.input:
        input_csv = Path(args.input)
    else:
        input_csv = find_latest_run_csv(output_dir)
        if not input_csv:
            print("未找到可导出的 CSV，请先运行 nnaa_fetch.py run 或 ai", file=sys.stderr)
            return 1

    stats = export_key_results(input_csv, output_dir, run_label=args.label)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


def cmd_run(args) -> int:
    """完整流水线：search → ai → export"""
    config_path = Path(args.config)
    cfg = _load_cfg(config_path)
    ai_cfg = cfg.get("ai_screening") or {}
    ai_enabled = bool(ai_cfg.get("enabled", False)) or args.force_ai

    if args.dry_run:
        return cmd_search(args)

    print("\n>>> 阶段 1/3：文献检索与去重")
    if cmd_search(args) != 0:
        return 1

    print("\n>>> 阶段 2/3：AI 文献判断")
    if not ai_enabled:
        print("提示: config.yaml 中 ai_screening.enabled=false，跳过 AI。")
        print("      启用 AI 后重新运行，或使用: nnaa_fetch.py run --force-ai")
    else:
        ai_tmp = None
        if args.force_ai and not ai_cfg.get("enabled", False):
            ai_tmp = _write_temp_cfg(cfg, {"ai_screening.enabled": True})
            ai_args = argparse.Namespace(
                config=str(ai_tmp),
                resume=args.resume,
                clear_checkpoint=args.clear_checkpoint,
                verbose=args.verbose,
            )
        else:
            ai_args = argparse.Namespace(
                config=args.config,
                resume=args.resume,
                clear_checkpoint=args.clear_checkpoint,
                verbose=args.verbose,
            )
        try:
            if cmd_ai(ai_args) != 0:
                return 1
        finally:
            if ai_tmp:
                ai_tmp.unlink(missing_ok=True)

    print("\n>>> 阶段 3/3：关键结果导出")
    export_args = argparse.Namespace(
        input=args.input,
        output_dir=args.output_dir,
        label=args.label,
    )
    if ai_enabled:
        return cmd_export(export_args)

    print("未启用 AI，跳过关键结果导出（关键文件需 ai_pass 标记）。")
    print("原始检索结果见 output/ 下日期目录。")
    return 0


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-c", "--config", default=str(DEFAULT_CONFIG), help="配置文件路径")
    common.add_argument("--date-from", metavar="YYYY-MM-DD", help="检索起始日期")
    common.add_argument("--date-to", metavar="YYYY-MM-DD", help="检索截止日期")
    common.add_argument(
        "--track",
        choices=["pathway", "enzymatic", "fermentation", "chemical", "hybrid", "all"],
        default="all",
        help="检索轨道",
    )
    common.add_argument(
        "--source",
        choices=["all", "pubmed", "crossref"],
        default="all",
        help="数据源（search/run 子命令）",
    )
    common.add_argument("-v", "--verbose", action="store_true")
    common.add_argument("--dry-run", action="store_true", help="仅预览检索式，不请求 API")
    common.add_argument("-i", "--input", help="export 时指定输入 CSV")
    common.add_argument("-o", "--output-dir", default=str(SCRIPT_DIR / "output"))
    common.add_argument("--label", help="关键结果子目录名")
    common.add_argument("--resume", action="store_true", help="AI 断点续跑")
    common.add_argument("--clear-checkpoint", action="store_true", help="清除 AI checkpoint")
    common.add_argument("--force-ai", action="store_true", help="run 时强制跑 AI（即使 config 未启用）")
    common.add_argument("--all-years", action="store_true", help="import-master: 不限日期，扫描全库")
    common.add_argument("--skip-ai", action="store_true", help="import-master: 跳过 AI 筛选")

    parser = argparse.ArgumentParser(
        description="非天然氨基酸文献流水线：检索 → AI 判断 → 关键结果输出",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
        parents=[common],
    )
    sub = parser.add_subparsers(dest="command", required=False)
    sub.add_parser("run", parents=[common], help="完整流水线：search → ai → export").set_defaults(func=cmd_run)
    sub.add_parser("search", parents=[common], help="仅检索与去重").set_defaults(func=cmd_search)
    sub.add_parser("ai", parents=[common], help="仅 AI 判断").set_defaults(func=cmd_ai)
    sub.add_parser("import-master", parents=[common], help="从酶库主库检索 NNAA 相关文献并并入 output").set_defaults(func=cmd_import_master)
    sub.add_parser("export", parents=[common], help="仅导出关键结果").set_defaults(func=cmd_export)
    sub.add_parser("keywords", help="打印检索关键词与检索式").set_defaults(func=cmd_keywords)
    parser.set_defaults(func=cmd_run)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

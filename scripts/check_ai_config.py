#!/usr/bin/env python3
"""检查 AI 筛选配置是否就绪（不调用 API）。"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def main():
    config_path = ROOT / "config.yaml"
    if not config_path.exists():
        print(f"缺少 {config_path}，请从 config.yaml.example 复制。", file=sys.stderr)
        sys.exit(1)

    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    ai = cfg.get("ai_screening") or {}
    key_env = ai.get("api_key_env", "NNAA_FETCH_OPENAI_API_KEY")
    has_key = bool(os.environ.get(key_env, "").strip() or (ai.get("api_key") or "").strip())
    api_base = (
        os.environ.get("NNAA_FETCH_OPENAI_BASE", "").strip()
        or ai.get("api_base", "https://api.openai.com/v1")
    )

    print("AI 筛选配置检查")
    print("-" * 40)
    print(f"enabled:              {ai.get('enabled', False)}")
    print(f"api_base:             {api_base}")
    print(f"model:                {ai.get('model', 'gpt-4o-mini')}")
    print(f"batch_size:           {ai.get('batch_size', 6)}")
    print(f"only_keep_passed:     {ai.get('only_keep_passed', True)}")
    print(f"keep_on_batch_error:  {ai.get('keep_on_batch_error', True)}")
    print(f"API Key ({key_env}): {'已配置' if has_key else '未配置'}")
    print("-" * 40)

    if not has_key:
        print("\n下一步：")
        print(f"  export {key_env}=sk-...")
        print("  或在 config.yaml 的 ai_screening.api_key 填写")
        print("\n试跑：")
        print("  python3 scripts/test_ai_screening.py -i output/2026-06-25/all.csv --limit 6 --keep-all")
        sys.exit(1)

    if not ai.get("enabled", False):
        print("\nKey 已就绪。正式启用请在 config.yaml 设 ai_screening.enabled: true")
    else:
        print("\n已启用。完整流水线：python3 nnaa_fetch.py run")

    print("\n试跑（不丢弃未通过项）：")
    print("  python3 scripts/test_ai_screening.py -i output/2026-06-25/all.csv --limit 6 --keep-all")


if __name__ == "__main__":
    main()

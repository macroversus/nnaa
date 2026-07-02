"""统一路径引导：确保项目根目录与 lib/ 均在 sys.path 中。"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
LIB_DIR = ROOT_DIR / "lib"

for path in (str(ROOT_DIR), str(LIB_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

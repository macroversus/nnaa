from __future__ import annotations

import hashlib
import json
import logging
from datetime import date, datetime
from pathlib import Path
from typing import List, Set

logger = logging.getLogger(__name__)

CHECKPOINT_FILE = "crossref_checkpoint.json"


def _signature(queries: List[str], date_from: date, date_to: date) -> str:
    raw = f"{date_from.isoformat()}|{date_to.isoformat()}|" + "\n".join(queries)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def load_checkpoint(state_dir: str, queries: List[str], date_from: date, date_to: date) -> Set[str]:
    path = Path(state_dir) / CHECKPOINT_FILE
    sig = _signature(queries, date_from, date_to)
    if not path.exists():
        return set()
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return set()
    if data.get("signature") != sig:
        logger.info("CrossRef checkpoint 与当前任务不匹配，将重新开始")
        return set()
    done = set(data.get("completed_queries") or [])
    if done:
        logger.info(f"CrossRef checkpoint: 已完成 {len(done)}/{len(queries)} 条检索式")
    return done


def save_checkpoint(
    state_dir: str,
    queries: List[str],
    date_from: date,
    date_to: date,
    completed: Set[str],
) -> None:
    Path(state_dir).mkdir(parents=True, exist_ok=True)
    path = Path(state_dir) / CHECKPOINT_FILE
    payload = {
        "signature": _signature(queries, date_from, date_to),
        "date_from": date_from.isoformat(),
        "date_to": date_to.isoformat(),
        "total_queries": len(queries),
        "completed_queries": sorted(completed),
        "updated_at": datetime.now().isoformat(timespec="seconds"),
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def clear_checkpoint(state_dir: str) -> None:
    path = Path(state_dir) / CHECKPOINT_FILE
    if path.exists():
        path.unlink()
        logger.info(f"已清除 CrossRef checkpoint: {path}")

from __future__ import annotations

import hashlib
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from deduplicator import normalize_doi, sanitize_record_for_csv

logger = logging.getLogger(__name__)

PENDING_FILE = "pending_ai.json"
CHECKPOINT_FILE = "ai_checkpoint.json"


def records_signature(records: List[Dict]) -> str:
    dois = sorted(
        {normalize_doi(r.get("doi") or "") for r in records if normalize_doi(r.get("doi") or "")}
    )
    raw = "\n".join(dois)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def save_pending(
    state_dir: str,
    records: List[Dict],
    meta: Optional[Dict[str, Any]] = None,
) -> str:
    Path(state_dir).mkdir(parents=True, exist_ok=True)
    path = Path(state_dir) / PENDING_FILE
    clean_records = [sanitize_record_for_csv(dict(r)) for r in records]
    payload = {
        "signature": records_signature(clean_records),
        "record_count": len(clean_records),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "meta": meta or {},
        "records": clean_records,
    }
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    logger.info(f"去重后待 AI 列表已保存: {path} ({len(clean_records)} 条)")
    return str(path)


def load_pending(state_dir: str) -> Tuple[List[Dict], Dict[str, Any]]:
    path = Path(state_dir) / PENDING_FILE
    if not path.exists():
        raise FileNotFoundError(f"待 AI 列表不存在: {path}（请先完整跑抓取+去重，或检查 state_dir）")
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    records = data.get("records") or []
    if not records:
        raise ValueError(f"待 AI 列表为空: {path}")
    meta = {
        "signature": data.get("signature", ""),
        "record_count": data.get("record_count", len(records)),
        "created_at": data.get("created_at", ""),
        **(data.get("meta") or {}),
    }
    sig = records_signature(records)
    if meta.get("signature") and meta["signature"] != sig:
        meta["signature"] = sig
        logger.warning("pending_ai.json 内 signature 与记录不一致，已按当前记录重算")
    else:
        meta["signature"] = sig
    logger.info(f"加载待 AI 列表: {path} ({len(records)} 条)")
    return records, meta


def pending_exists(state_dir: str) -> bool:
    return (Path(state_dir) / PENDING_FILE).exists()


def load_ai_checkpoint(
    state_dir: str,
    expected_signature: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    path = Path(state_dir) / CHECKPOINT_FILE
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        logger.warning("AI checkpoint 读取失败，将重新开始")
        return None
    if expected_signature and data.get("signature") != expected_signature:
        logger.info("AI checkpoint 与当前待筛列表不匹配，将重新开始")
        return None
    results = data.get("results_by_doi") or {}
    failed = data.get("batch_failed_dois") or []
    done = len(results) + len(failed)
    total = data.get("total_records") or 0
    logger.info(f"AI checkpoint: 已处理 {done}/{total} 条")
    return data


def save_ai_checkpoint(
    state_dir: str,
    signature: str,
    total_records: int,
    batch_size: int,
    results_by_doi: Dict[str, Dict],
    batch_failed_dois: Set[str],
    completed_batches: int,
) -> None:
    Path(state_dir).mkdir(parents=True, exist_ok=True)
    path = Path(state_dir) / CHECKPOINT_FILE
    payload = {
        "signature": signature,
        "total_records": total_records,
        "batch_size": batch_size,
        "completed_batches": completed_batches,
        "results_by_doi": results_by_doi,
        "batch_failed_dois": sorted(batch_failed_dois),
        "updated_at": datetime.now().isoformat(timespec="seconds"),
    }
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def clear_ai_checkpoint(state_dir: str) -> None:
    path = Path(state_dir) / CHECKPOINT_FILE
    if path.exists():
        path.unlink()
        logger.info(f"已清除 AI checkpoint: {path}")


def clear_pending(state_dir: str) -> None:
    path = Path(state_dir) / PENDING_FILE
    if path.exists():
        path.unlink()
        logger.info(f"已清除待 AI 列表: {path}")

from __future__ import annotations

import logging
from typing import List, Set

from deduplicator import load_existing_dois, normalize_doi

logger = logging.getLogger(__name__)


def load_local_pdf_dois(csv_path: str, doi_column: str = "DOI") -> Set[str]:
    if not csv_path or not str(csv_path).strip():
        return set()
    return load_existing_dois(csv_path, doi_column)


def enrich_records(records: List[dict], local_dois: Set[str]) -> None:
    for rec in records:
        doi = normalize_doi(rec.get("doi") or "")
        rec["local_has_pdf"] = "true" if doi in local_dois else "false"

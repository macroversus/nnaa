"""CrossRef 检索式构建 — 关键词来自 nnaa_keywords.py"""

from __future__ import annotations

from typing import List, Optional

from nnaa_keywords import build_crossref_queries

TRACK_QUERY_BUILDERS = {
    "pathway": lambda max_queries=None: build_crossref_queries("pathway", max_queries),
    "enzymatic": lambda max_queries=None: build_crossref_queries("enzymatic", max_queries),
    "fermentation": lambda max_queries=None: build_crossref_queries("fermentation", max_queries),
    "chemical": lambda max_queries=None: build_crossref_queries("chemical", max_queries),
    "hybrid": lambda max_queries=None: build_crossref_queries("hybrid", max_queries),
    "gce": lambda max_queries=None: build_crossref_queries("gce", max_queries),
}


def build_pathway_queries(max_queries: Optional[int] = None) -> List[str]:
    return build_crossref_queries("pathway", max_queries)


def build_enzymatic_queries(max_queries: Optional[int] = None) -> List[str]:
    return build_crossref_queries("enzymatic", max_queries)


def build_fermentation_queries(max_queries: Optional[int] = None) -> List[str]:
    return build_crossref_queries("fermentation", max_queries)


def build_chemical_queries(max_queries: Optional[int] = None) -> List[str]:
    return build_crossref_queries("chemical", max_queries)


def build_hybrid_queries(max_queries: Optional[int] = None) -> List[str]:
    return build_crossref_queries("hybrid", max_queries)

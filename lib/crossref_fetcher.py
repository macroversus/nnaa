import time
import logging
import requests
from datetime import datetime, date
from typing import List, Dict, Optional, Set

logger = logging.getLogger(__name__)

CROSSREF_WORKS_URL = "https://api.crossref.org/works"


class CrossRefFetcher:
    def __init__(self, config: dict):
        self.cfg = config.get("crossref", {})
        self.mailto = self.cfg.get("mailto", "")
        self.rows_per_query = int(self.cfg.get("rows_per_query", 200))
        self.filter_types = self.cfg.get("filter_types", ["journal-article"])
        _qp = (self.cfg.get("query_param") or "query").strip().lower()
        if _qp in ("bibliographic", "query.bibliographic"):
            self._query_api_key = "query.bibliographic"
        elif _qp in ("title", "query.title"):
            self._query_api_key = "query.title"
        else:
            self._query_api_key = "query"
        default_cap = self.rows_per_query * 20
        self.limit_per_query = int(self.cfg.get("limit_per_query", default_cap))
        if self.cfg.get("request_interval_sec") is not None:
            self.request_interval = float(self.cfg["request_interval_sec"])
        else:
            self.request_interval = 0.2 if self.mailto else 1.0

    def _base_headers(self) -> dict:
        headers = {"User-Agent": "LiteratureFetcher/1.0"}
        if self.mailto:
            headers["User-Agent"] += f" (mailto:{self.mailto})"
        return headers

    def _get(self, url: str, params: dict, retries: int = 3) -> Optional[requests.Response]:
        for attempt in range(retries):
            try:
                resp = requests.get(url, params=params,
                                    headers=self._base_headers(), timeout=30)
                if resp.status_code == 429:
                    wait = int(resp.headers.get("Retry-After", 60))
                    logger.warning(f"CrossRef 429 限流，等待 {wait}s 后重试")
                    time.sleep(wait)
                    continue
                if resp.status_code in (404, 410):
                    return resp
                resp.raise_for_status()
                time.sleep(self.request_interval)
                return resp
            except requests.RequestException as e:
                wait = 2 ** attempt
                logger.warning(f"CrossRef 请求失败 (attempt {attempt+1}/{retries}): {e}，{wait}s 后重试")
                time.sleep(wait)
        return None

    def _build_filter(self, date_from: date, date_to: date) -> str:
        filters = [
            f"from-pub-date:{date_from.strftime('%Y-%m-%d')}",
            f"until-pub-date:{date_to.strftime('%Y-%m-%d')}",
        ]
        if self.filter_types:
            for t in self.filter_types:
                filters.append(f"type:{t}")
        return ",".join(filters)

    def search(self, query: str, date_from: date, date_to: date) -> List[Dict]:
        filter_str = self._build_filter(date_from, date_to)

        params = {
            self._query_api_key: query,
            "filter": filter_str,
            "rows": 0,
            "select": "DOI",
        }
        resp = self._get(CROSSREF_WORKS_URL, params)
        if not resp:
            logger.error(f"CrossRef 查询失败: {query}")
            return []

        data = resp.json()
        total = data.get("message", {}).get("total-results", 0)
        if total == 0:
            logger.info(f"CrossRef 查询无结果: {query}")
            return []

        logger.info(f"CrossRef 查询 '{query[:60]}' 共 {total} 条，日期 {date_from}~{date_to}")

        results = []
        cursor = "*"
        fetched = 0
        max_fetch = min(total, self.limit_per_query)

        while fetched < max_fetch:
            batch_size = min(self.rows_per_query, max_fetch - fetched)
            params = {
                self._query_api_key: query,
                "filter": filter_str,
                "rows": batch_size,
                "cursor": cursor,
                "select": "DOI,title,container-title,published,type",
            }
            resp = self._get(CROSSREF_WORKS_URL, params)
            if not resp:
                logger.warning(f"CrossRef 分页失败，停止该查询 cursor={cursor}")
                break

            msg = resp.json().get("message", {})
            items = msg.get("items", [])
            if not items:
                break

            for item in items:
                record = self._parse_item(item)
                if record["doi"]:
                    results.append(record)

            fetched += len(items)
            cursor = msg.get("next-cursor", "")
            if not cursor:
                break

            logger.debug(f"  CrossRef 已获取 {fetched} 条")

        return results

    @staticmethod
    def _unescape(text: str) -> str:
        """解码 CrossRef 返回的 HTML 实体（&amp; → &, &lt; → < 等）。"""
        import html
        return html.unescape(text) if text else text

    def _parse_item(self, item: dict) -> Dict:
        doi = (item.get("DOI") or "").strip()

        titles = item.get("title", [])
        title = self._unescape(titles[0]) if titles else ""

        containers = item.get("container-title", [])
        journal = self._unescape(containers[0]) if containers else ""

        pub_date = ""
        pub = item.get("published") or item.get("published-print") or item.get("published-online")
        if pub:
            parts = pub.get("date-parts", [[]])[0]
            pub_date = "-".join(str(p) for p in parts if p)

        return {
            "source": "crossref",
            "pmid": "",
            "doi": doi,
            "title": title,
            "journal": journal,
            "pub_date": pub_date,
            "abstract": "",
            "publication_types": (item.get("type") or "").replace("-", " "),
        }

    def fetch_since(
        self,
        queries: List[str],
        date_from: date,
        date_to: Optional[date] = None,
        *,
        state_dir: Optional[str] = None,
        resume: bool = False,
    ) -> List[Dict]:
        if not self.cfg.get("enabled", True):
            logger.info("CrossRef 数据源已禁用，跳过")
            return []

        date_to = date_to or date.today()
        all_results: List[Dict] = []
        total_q = len(queries)
        completed: Set[str] = set()

        if state_dir and resume:
            from crossref_checkpoint import load_checkpoint, save_checkpoint

            completed = load_checkpoint(state_dir, queries, date_from, date_to)
        elif state_dir and not resume:
            from crossref_checkpoint import clear_checkpoint

            clear_checkpoint(state_dir)

        pending = [q for q in queries if q not in completed]
        if completed:
            logger.info(f"CrossRef 跳过已完成检索式 {len(completed)} 条，剩余 {len(pending)} 条")

        for idx, query in enumerate(pending, 1):
            done_count = len(completed) + idx
            if total_q >= 10 and (idx == 1 or idx % 10 == 0 or idx == len(pending)):
                logger.info(f"CrossRef 进度: {done_count}/{total_q}")
            records = self.search(query, date_from, date_to)
            all_results.extend(records)
            if state_dir:
                from crossref_checkpoint import save_checkpoint

                completed.add(query)
                save_checkpoint(state_dir, queries, date_from, date_to, completed)

        if state_dir and len(completed) >= total_q:
            from crossref_checkpoint import clear_checkpoint

            logger.info("CrossRef 全部检索式已完成，清除 checkpoint")
            clear_checkpoint(state_dir)

        seen: set = set()
        unique: List[Dict] = []
        for r in all_results:
            doi_key = r["doi"].lower()
            if doi_key not in seen:
                seen.add(doi_key)
                unique.append(r)

        logger.info(f"CrossRef 全部查询合并后共 {len(unique)} 个唯一 DOI")
        return unique

    def fetch_work_by_doi(self, doi: str) -> Optional[Dict]:
        doi = (doi or "").strip()
        if not doi:
            return None
        url = f"{CROSSREF_WORKS_URL}/{doi}"
        resp = self._get(url, {})
        if not resp or resp.status_code in (404, 410):
            return None
        item = resp.json().get("message") or {}
        if not item.get("DOI"):
            return None
        return self._parse_item(item)

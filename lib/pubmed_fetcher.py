import time
import logging
import requests
import xml.etree.ElementTree as ET
from datetime import datetime, date
from typing import List, Dict, Optional

logger = logging.getLogger(__name__)

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class PubMedFetcher:
    def __init__(self, config: dict):
        self.cfg = config.get("pubmed", {})
        self.email   = self.cfg.get("email", "")
        self.api_key = self.cfg.get("api_key", "")
        self.batch_size = int(self.cfg.get("batch_size", 200))
        self.efetch_batch_size = int(
            self.cfg.get("efetch_batch_size") or min(self.batch_size, 200)
        )
        self.max_pmids = int(self.cfg.get("max_pmids") or 0)
        self.request_interval = 0.11 if self.api_key else 0.34

    def _base_params(self) -> dict:
        params = {"tool": "lit_fetcher", "email": self.email, "retmode": "json"}
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def _get(self, url: str, params: dict, retries: int = 5) -> Optional[requests.Response]:
        for attempt in range(retries):
            try:
                resp = requests.get(url, params=params, timeout=30)
                resp.raise_for_status()
                time.sleep(self.request_interval)
                return resp
            except requests.RequestException as e:
                wait = 2 ** attempt
                logger.warning(f"PubMed 请求失败 (attempt {attempt+1}/{retries}): {e}，{wait}s 后重试")
                time.sleep(wait)
        return None

    def _get_json(self, url: str, params: dict, retries: int = 5) -> Optional[dict]:
        """带 JSON 解析重试的 GET：应对 PubMed 偶发返回带控制字符的无效响应。"""
        for attempt in range(retries):
            try:
                resp = requests.get(url, params=params, timeout=30)
                resp.raise_for_status()
                time.sleep(self.request_interval)
                return resp.json()
            except (requests.RequestException, requests.exceptions.JSONDecodeError,
                    ValueError) as e:
                wait = 2 ** attempt
                logger.warning(f"PubMed 请求/解析失败 (attempt {attempt+1}/{retries}): {e}，{wait}s 后重试")
                time.sleep(wait)
        logger.error(f"PubMed 请求在 {retries} 次重试后仍失败，跳过")
        return None

    def _post(self, url: str, data: dict, retries: int = 3) -> Optional[requests.Response]:
        for attempt in range(retries):
            try:
                resp = requests.post(url, data=data, timeout=60)
                resp.raise_for_status()
                time.sleep(self.request_interval)
                return resp
            except requests.RequestException as e:
                wait = 2 ** attempt
                logger.warning(f"PubMed POST 请求失败 (attempt {attempt+1}/{retries}): {e}，{wait}s 后重试")
                time.sleep(wait)
        return None

    def search_pmids(self, query: str, date_from: date, date_to: date) -> List[str]:
        """搜索 PMID 列表。
        当结果数超过 PubMed esearch 的硬限制（9999）时，自动按年拆分递归检索。
        """
        return self._search_pmids_range(query, date_from, date_to)

    def _search_pmids_range(self, query: str, date_from: date, date_to: date) -> List[str]:
        """单一日期段检索；超过 9999 条时拆年递归。"""
        mindate = date_from.strftime("%Y/%m/%d")
        maxdate = date_to.strftime("%Y/%m/%d")

        # 先只取 count，不下载 ID
        params = {
            **self._base_params(),
            "db": "pubmed",
            "term": query,
            "mindate": mindate,
            "maxdate": maxdate,
            "datetype": "pdat",
            "retmax": 0,
            "retstart": 0,
            "retmode": "json",
        }
        data = self._get_json(ESEARCH_URL, params)
        if not data:
            logger.warning(f"无法获取 PubMed 结果数: {query[:50]}...")
            return []

        total = int(data.get("esearchresult", {}).get("count", 0))
        if total == 0:
            logger.info(f"PubMed 查询无结果: {query[:80]}")
            return []

        # PubMed esearch 硬限制：retstart 不能超过 9999
        ESEARCH_MAX = 9999
        if total > ESEARCH_MAX:
            # 按年拆分；若同一年还超限则按月拆
            return self._split_and_search(query, date_from, date_to, total)

        logger.info(f"PubMed 查询 '{query[:60]}' 共 {total} 条，日期 {mindate}~{maxdate}")
        if self.max_pmids and total > self.max_pmids:
            logger.warning(f"  结果数 {total} 超过 max_pmids={self.max_pmids}，截断")
            total = self.max_pmids

        # 分页下载 PMID
        all_pmids: List[str] = []
        for start in range(0, total, self.batch_size):
            params.update({"retmax": min(self.batch_size, total - start), "retstart": start})
            page_data = self._get_json(ESEARCH_URL, params)
            if not page_data:
                logger.warning(f"分页获取失败，跳过 start={start}")
                continue
            ids = page_data.get("esearchresult", {}).get("idlist", [])
            all_pmids.extend(ids)
            logger.debug(f"  已获取 {len(all_pmids)}/{total}")

        return all_pmids

    def _split_and_search(self, query: str, date_from: date, date_to: date,
                          total: int) -> List[str]:
        """将日期范围按年（或按月）拆分递归检索，绕过 esearch 9999 上限。"""
        from datetime import timedelta
        import calendar

        span_days = (date_to - date_from).days
        all_pmids: List[str] = []

        if span_days >= 365:
            # 按年拆分
            year = date_from.year
            end_year = date_to.year
            while year <= end_year:
                seg_from = date(year, 1, 1) if year > date_from.year else date_from
                seg_to   = date(year, 12, 31) if year < date_to.year else date_to
                seg_pmids = self._search_pmids_range(query, seg_from, seg_to)
                all_pmids.extend(seg_pmids)
                year += 1
        elif span_days >= 28:
            # 按月拆分
            cur = date_from
            while cur <= date_to:
                last_day = calendar.monthrange(cur.year, cur.month)[1]
                seg_from = cur
                seg_to   = min(date(cur.year, cur.month, last_day), date_to)
                seg_pmids = self._search_pmids_range(query, seg_from, seg_to)
                all_pmids.extend(seg_pmids)
                # 下一个月第一天
                if cur.month == 12:
                    cur = date(cur.year + 1, 1, 1)
                else:
                    cur = date(cur.year, cur.month + 1, 1)
        else:
            # 范围已经很小但还超限，截断处理
            logger.warning(f"  日期段 {date_from}~{date_to} 结果仍超 9999，截断至 9999")
            params = {
                **self._base_params(),
                "db": "pubmed",
                "term": query,
                "mindate": date_from.strftime("%Y/%m/%d"),
                "maxdate": date_to.strftime("%Y/%m/%d"),
                "datetype": "pdat",
                "retmax": 9999,
                "retstart": 0,
                "retmode": "json",
            }
            page_data = self._get_json(ESEARCH_URL, params)
            if page_data:
                all_pmids.extend(page_data.get("esearchresult", {}).get("idlist", []))

        return all_pmids

    def _efetch_batch(self, batch: List[str]) -> List[Dict]:
        data = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "xml",
            "retmode": "xml",
            "tool": "lit_fetcher",
            "email": self.email,
        }
        if self.api_key:
            data["api_key"] = self.api_key

        resp = self._post(EFETCH_URL, data)
        if not resp:
            if len(batch) <= 1:
                logger.warning(f"efetch 失败，跳过 PMID {batch[0] if batch else '?'}")
                return []
            mid = len(batch) // 2
            logger.warning(f"efetch 失败，拆分批次 {len(batch)} -> {mid}+{len(batch) - mid}")
            return self._efetch_batch(batch[:mid]) + self._efetch_batch(batch[mid:])

        return self._parse_pubmed_xml(resp.text)

    def fetch_details(self, pmids: List[str]) -> List[Dict]:
        if not pmids:
            return []

        results = []
        total = len(pmids)
        for i in range(0, total, self.efetch_batch_size):
            batch = pmids[i: i + self.efetch_batch_size]
            parsed = self._efetch_batch(batch)
            results.extend(parsed)
            if (i + self.efetch_batch_size) % 500 == 0 or i + self.efetch_batch_size >= total:
                logger.info(f"  PubMed efetch 进度: {min(i + len(batch), total)}/{total}")

        return results

    def _parse_pubmed_xml(self, xml_text: str) -> List[Dict]:
        records = []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as e:
            logger.error(f"XML 解析错误: {e}")
            return records

        for article in root.findall(".//PubmedArticle"):
            record: Dict = {
                "source": "pubmed",
                "pmid": "",
                "doi": "",
                "title": "",
                "journal": "",
                "pub_date": "",
                "abstract": "",
                "publication_types": "",
            }

            pmid_el = article.find(".//PMID")
            if pmid_el is not None:
                record["pmid"] = pmid_el.text or ""

            for aid in article.findall(".//ArticleId"):
                if aid.get("IdType") == "doi":
                    record["doi"] = (aid.text or "").strip()
                    break

            title_el = article.find(".//ArticleTitle")
            if title_el is not None:
                record["title"] = "".join(title_el.itertext()).strip()

            journal_el = article.find(".//Journal/Title")
            if journal_el is not None:
                record["journal"] = journal_el.text or ""

            pt_nodes = article.findall(".//PublicationTypeList/PublicationType")
            if pt_nodes:
                record["publication_types"] = "; ".join(
                    "".join(n.itertext()).strip() for n in pt_nodes if n is not None
                )

            abs_parts = []
            for at in article.findall(".//Abstract/AbstractText"):
                label = at.get("Label") or ""
                chunk = "".join(at.itertext()).strip()
                if not chunk:
                    continue
                if label:
                    abs_parts.append(f"{label}: {chunk}")
                else:
                    abs_parts.append(chunk)
            if abs_parts:
                record["abstract"] = "\n".join(abs_parts)

            pub_date = self._extract_date(article)
            record["pub_date"] = pub_date

            records.append(record)

        return records

    def _extract_date(self, article_el) -> str:
        for date_path in [".//PubDate", ".//ArticleDate"]:
            d = article_el.find(date_path)
            if d is not None:
                year  = (d.findtext("Year")  or "").strip()
                month = (d.findtext("Month") or "").strip()
                day   = (d.findtext("Day")   or "").strip()
                if year:
                    parts = [year, month, day]
                    return "-".join(p for p in parts if p)
        return ""

    def fetch_since(self, queries: List[str], date_from: date, date_to: Optional[date] = None) -> List[Dict]:
        if not self.cfg.get("enabled", True):
            logger.info("PubMed 数据源已禁用，跳过")
            return []

        date_to = date_to or date.today()
        all_pmids: List[str] = []

        for query in queries:
            pmids = self.search_pmids(query, date_from, date_to)
            all_pmids.extend(pmids)

        unique_pmids = list(dict.fromkeys(all_pmids))
        logger.info(f"PubMed 全部查询合并后共 {len(unique_pmids)} 个唯一 PMID")

        if not unique_pmids:
            return []

        if self.max_pmids > 0 and len(unique_pmids) > self.max_pmids:
            logger.warning(
                f"PubMed PMID 数 {len(unique_pmids)} 超过 max_pmids={self.max_pmids}，截断"
            )
            unique_pmids = unique_pmids[: self.max_pmids]

        return self.fetch_details(unique_pmids)

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

    def _get(self, url: str, params: dict, retries: int = 3) -> Optional[requests.Response]:
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
        mindate = date_from.strftime("%Y/%m/%d")
        maxdate = date_to.strftime("%Y/%m/%d")

        params = {
            **self._base_params(),
            "db": "pubmed",
            "term": query,
            "mindate": mindate,
            "maxdate": maxdate,
            "datetype": "pdat",
            "retmax": 0,
            "retstart": 0,
        }
        params.pop("retmode", None)
        params["retmode"] = "json"

        resp = self._get(ESEARCH_URL, params)
        if not resp:
            logger.error(f"无法获取 PubMed 搜索结果: {query}")
            return []

        data = resp.json()
        total = int(data.get("esearchresult", {}).get("count", 0))
        if total == 0:
            logger.info(f"PubMed 查询无结果: {query}")
            return []

        logger.info(f"PubMed 查询 '{query[:60]}' 共 {total} 条，日期 {mindate}~{maxdate}")

        all_pmids: List[str] = []
        for start in range(0, total, self.batch_size):
            params.update({"retmax": self.batch_size, "retstart": start})
            resp = self._get(ESEARCH_URL, params)
            if not resp:
                logger.warning(f"分页获取失败，跳过 start={start}")
                continue
            ids = resp.json().get("esearchresult", {}).get("idlist", [])
            all_pmids.extend(ids)
            logger.debug(f"  已获取 {len(all_pmids)}/{total}")

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

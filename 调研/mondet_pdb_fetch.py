#!/usr/bin/env python3
"""
MONDE·T PDB → 论文反查脚本

流程：
1. 读取 MONDE·T CSV，提取所有 PDB 结构 ID
2. 调用 RCSB PDB REST API 查询每个结构对应的 PMID/DOI
3. 用 PubMed efetch 批量获取标题/摘要
4. 与现有 all.csv 对比，找出新增文献
5. 输出到 /root/非天然氨基酸/调研/mondet_new_papers.csv

用法：python3 mondet_pdb_fetch.py [--limit N] [--resume]
"""
import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

BASE_DIR = Path("/root/非天然氨基酸/调研")
ALL_CSV = Path("/root/非天然氨基酸/output/2026-07-09_190703/all.csv")
OUT_CSV = BASE_DIR / "mondet_new_papers.csv"
CHECKPOINT = BASE_DIR / "mondet_pdb_checkpoint.json"

RCSB_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
PM_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PM_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def load_existing_dois():
    """Load DOIs already in our all.csv"""
    dois = set()
    if ALL_CSV.exists():
        with open(ALL_CSV, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                doi_raw = (row.get('doi') or '').strip().lower()
                doi = re.sub(r'^https?://(dx\.)?doi\.org/', '', doi_raw)
                if doi:
                    dois.add(doi)
    print(f"Existing DOIs in all.csv: {len(dois)}")
    return dois


def load_mondet_pdb_ids():
    """Extract all unique PDB IDs from MONDE·T"""
    pdb_ids = set()
    mondet_file = BASE_DIR / "mondet_tab.csv"
    with open(mondet_file, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entity_ids = row.get('entity_IDs', '')
            ids = re.findall(r'[A-Z0-9]{4}', entity_ids)
            for eid in ids:
                if len(eid) == 4:
                    pdb_ids.add(eid)
    print(f"Unique PDB IDs from MONDE·T: {len(pdb_ids)}")
    return pdb_ids


def fetch_pdb_pmid(pdb_id: str) -> dict:
    """Fetch PMID and DOI for a PDB entry via RCSB REST API"""
    url = RCSB_URL.format(pdb_id=pdb_id)
    req = urllib.request.Request(url, headers={"User-Agent": "NNAA-research/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            data = json.load(r)
        result = {'pdb_id': pdb_id, 'pmid': '', 'doi': '', 'title': '', 'year': ''}
        # Primary citation
        citation = data.get('rcsb_primary_citation') or {}
        pmid = str(citation.get('pdbx_database_id_PubMed') or '').strip()
        doi = str(citation.get('pdbx_database_id_DOI') or '').strip()
        result['pmid'] = pmid
        result['doi'] = doi.lower()
        result['title'] = str(citation.get('title') or '').strip()
        result['year'] = str(citation.get('year') or '').strip()
        return result
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {'pdb_id': pdb_id, 'error': '404'}
        raise
    except Exception as e:
        return {'pdb_id': pdb_id, 'error': str(e)[:80]}


def fetch_pubmed_abstracts(pmids: list) -> dict:
    """Batch fetch abstracts from PubMed for a list of PMIDs"""
    if not pmids:
        return {}
    results = {}
    # Fetch in batches of 100
    for i in range(0, len(pmids), 100):
        batch = pmids[i:i+100]
        params = urllib.parse.urlencode({
            'db': 'pubmed',
            'id': ','.join(batch),
            'retmode': 'xml',
            'rettype': 'abstract',
        })
        try:
            import urllib.parse
            url = f"{PM_EFETCH}?{params}"
            req = urllib.request.Request(url, headers={"User-Agent": "NNAA-research/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                xml = r.read().decode('utf-8', errors='ignore')
            # Parse titles and abstracts from XML
            titles = re.findall(r'<ArticleTitle>([^<]+)</ArticleTitle>', xml)
            abstracts = re.findall(r'<AbstractText[^>]*>([^<]+)</AbstractText>', xml)
            pmid_matches = re.findall(r'<PMID[^>]*>(\d+)</PMID>', xml)
            for j, pmid in enumerate(pmid_matches[:len(titles)]):
                results[pmid] = {
                    'title': titles[j] if j < len(titles) else '',
                    'abstract': abstracts[j] if j < len(abstracts) else '',
                }
            time.sleep(0.35)
        except Exception as e:
            print(f"  PubMed efetch error for batch {i//100+1}: {e}")
    return results


def main():
    import urllib.parse
    parser = argparse.ArgumentParser(description="MONDE·T PDB → Paper lookup")
    parser.add_argument('--limit', type=int, default=0, help='Max PDB IDs to process (0=all)')
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    args = parser.parse_args()

    existing_dois = load_existing_dois()
    pdb_ids = load_mondet_pdb_ids()
    pdb_ids = sorted(pdb_ids)
    if args.limit:
        pdb_ids = pdb_ids[:args.limit]
        print(f"Limited to first {args.limit} PDB IDs")

    # Load checkpoint
    done_pdb = {}
    if args.resume and CHECKPOINT.exists():
        with open(CHECKPOINT) as f:
            done_pdb = json.load(f)
        print(f"Resumed from checkpoint: {len(done_pdb)} PDB IDs already processed")

    # Process PDB IDs
    new_papers = {}  # doi/pmid -> paper info
    errors = 0
    for i, pdb_id in enumerate(pdb_ids):
        if pdb_id in done_pdb:
            res = done_pdb[pdb_id]
        else:
            res = fetch_pdb_pmid(pdb_id)
            done_pdb[pdb_id] = res
            if i % 200 == 0:
                with open(CHECKPOINT, 'w') as f:
                    json.dump(done_pdb, f)
                print(f"  [{i}/{len(pdb_ids)}] checkpoint saved, errors={errors}")
            time.sleep(0.05)

        if res.get('error'):
            errors += 1
            continue

        doi = res.get('doi', '').lower()
        doi = re.sub(r'^https?://(dx\.)?doi\.org/', '', doi)
        pmid = res.get('pmid', '')

        # Check if already in our database
        if doi and doi not in existing_dois:
            key = doi or pmid
            if key and key not in new_papers:
                new_papers[key] = {
                    'pdb_id': pdb_id,
                    'doi': doi,
                    'pmid': pmid,
                    'title': res.get('title', ''),
                    'year': res.get('year', ''),
                    'source': 'mondet_pdb',
                }
        elif pmid and doi in existing_dois:
            pass  # already have it

    print(f"\nProcessed {len(pdb_ids)} PDB IDs")
    print(f"Errors: {errors}")
    print(f"New papers not in all.csv: {len(new_papers)}")

    # Save checkpoint final
    with open(CHECKPOINT, 'w') as f:
        json.dump(done_pdb, f)

    # Write output
    if new_papers:
        fieldnames = ['pdb_id', 'doi', 'pmid', 'title', 'year', 'source']
        with open(OUT_CSV, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for paper in new_papers.values():
                writer.writerow({k: paper.get(k, '') for k in fieldnames})
        print(f"Saved to {OUT_CSV}")

    # Summary stats
    years = {}
    for p in new_papers.values():
        y = p.get('year', '')
        if y:
            years[y] = years.get(y, 0) + 1
    print("Year distribution (top 10):")
    for y, c in sorted(years.items(), reverse=True)[:10]:
        print(f"  {y}: {c}")


if __name__ == '__main__':
    main()

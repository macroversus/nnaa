# 非天然氨基酸文献获取

独立项目，流程：**检索 → AI 判断 → 关键结果输出**

```bash
cd /root/非天然氨基酸
pip install -r requirements.txt
cp config.yaml.example config.yaml

python3 nnaa_fetch.py keywords       # 查看关键词
python3 nnaa_fetch.py run --dry-run  # 预览
python3 nnaa_fetch.py run            # 完整流水线
```

## 目录结构

```
非天然氨基酸/
├── nnaa_fetch.py           # ★ 唯一入口
├── nnaa_keywords.py        # ★ 检索关键词库
├── nnaa_queries.py           # 五轨道查询路由
├── nnaa_query_strategy.py    # CrossRef 检索式
├── export_key_results.py     # ★ 关键 CSV 导出
├── bootstrap.py              # 路径引导
├── config.yaml.example
├── lib/                      # 检索与 AI 管线
│   ├── fetch_literature.py   # 抓取编排
│   ├── ai_screener.py        # AI 筛选
│   ├── pubmed_fetcher.py
│   ├── crossref_fetcher.py
│   └── deduplicator.py       # …等
├── scripts/                  # 辅助工具
│   ├── check_ai_config.py
│   ├── test_ai_screening.py
│   └── repair_output_csv.py
├── docs/
│   ├── 使用说明.md
│   └── 架构说明.md
├── output/                   # 运行结果
│   └── key_results/          # ★ 关键交付 CSV
├── state/
└── logs/
```

## 文档

- [使用说明](docs/使用说明.md)
- [架构说明](docs/架构说明.md)

## 关键输出

`output/key_results/` 下的 `关键_*.csv` 为 AI 筛选通过后的最终交付物。

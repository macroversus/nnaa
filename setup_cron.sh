#!/bin/bash
# 每月 1 日 02:00 跑非天然氨基酸文献流水线
# 用法: bash setup_cron.sh

DIR="$(cd "$(dirname "$0")" && pwd)"
CRON_LINE="0 2 1 * * cd ${DIR} && /usr/bin/python3 nnaa_fetch.py run >> logs/cron.log 2>&1"

if crontab -l 2>/dev/null | grep -F "nnaa_fetch.py run" >/dev/null; then
  echo "cron 已存在 nnaa_fetch 任务"
else
  (crontab -l 2>/dev/null; echo "$CRON_LINE") | crontab -
  echo "已添加: $CRON_LINE"
fi

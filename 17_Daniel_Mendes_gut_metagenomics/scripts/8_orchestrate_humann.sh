#!/usr/bin/env bash
# 等 vJun23 db 下载安装完成 → 自动启动 HUMAnN → 等其跑完。
P=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
DB=/Work_bio/references/Metagenomics/metaphlan_vJun23
# 1) 等下载安装完成（mpavj23 session 结束 且 关键文件就位）
for i in $(seq 1 180); do   # 最多等 3h
  if ! tmux has-session -t mpavj23 2>/dev/null; then
    if ls "$DB"/${DB##*/} 2>/dev/null || ls "$DB"/mpa_vJun23*.pkl 2>/dev/null || ls "$DB"/mpa_vJun23*.bt2l 2>/dev/null; then
      echo "[$(date)] vJun23 db ready"; break
    fi
    echo "[$(date)] mpavj23 session ended but db files check..."; ls "$DB" | head; break
  fi
  sleep 60
done
# 2) 启动 HUMAnN
echo "[$(date)] launching HUMAnN"
tmux new-session -d -s humann17 "bash $P/scripts/7_run_humann.sh 2>&1 | tee $P/logs/humann_run.log"

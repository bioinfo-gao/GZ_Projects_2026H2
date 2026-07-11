#!/bin/bash
# ============================================================================
# Step 11 — 批量驱动:对全部 6 样跑完下游敲入验证分析(5→4→6→7,CD1A 另跑 8)
#   sarek(--tools tiddit)已于 2026-07-11 完成;本脚本把在 RAGH_153 试跑验证过的
#   下游定制分析推广到全 6 样。顺序按各脚本 header 依赖:5(拷贝数,产 mosdepth 基线)
#   必须最先,4/6/7 依赖它。串行逐样跑(每样 mosdepth -t6),避免与仍在跑的项目14
#   Mutect2 抢核。
#   用法: bash 11_run_all_downstream.sh [sample1 sample2 ...]  (缺省=全部6样)
# ============================================================================
set -uo pipefail
PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
cd "$PROJ/scripts"
SAMPLES=("$@")
[ ${#SAMPLES[@]} -eq 0 ] && SAMPLES=(CD1A_B125 RAGH_153 RAGH_273 MTTH_284 MTTH_412 MTTH_524)
LOG="$PROJ/logs/downstream_run_$(date +%m%d).log"
echo "==== downstream batch start $(date) | samples: ${SAMPLES[*]} ====" | tee -a "$LOG"

run_step(){  # step_script sample
  local scr="$1" smp="$2"
  echo ">>> [$(date +%H:%M:%S)] $scr $smp" | tee -a "$LOG"
  if bash "$scr" "$smp" >>"$LOG" 2>&1; then
    echo "    OK  $scr $smp" | tee -a "$LOG"
  else
    echo "    FAIL $scr $smp (见日志)" | tee -a "$LOG"; return 1
  fi
}

for S in "${SAMPLES[@]}"; do
  echo "======== SAMPLE $S ========" | tee -a "$LOG"
  run_step 5_copy_number.sh      "$S" || continue   # 必须最先(基线)
  run_step 4_integration_analysis.sh "$S"
  run_step 6_ki_integrity_check.sh   "$S"
  run_step 7_zygosity_analysis.sh    "$S"
  [ "$S" = "CD1A_B125" ] && run_step 8_cd1a_neo_status.sh "$S"
done
echo "==== downstream batch done $(date) ====" | tee -a "$LOG"

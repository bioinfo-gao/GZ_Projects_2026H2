#!/usr/bin/env bash
# Phase 1b — functional profiling (HUMAnN v3.9), assembly-free。
# 输入 = taxprofiler 去宿主后的 unmapped reads（work_taxprofiler 里，R1+R2 合并成单端喂 HUMAnN）。
# DB: ChocoPhlAn(16G) + UniRef90 full(34G) + MetaPhlAn SGB(vJan25)，humann_config 已就位。
# 并发：3 样本并行 × 18 threads = 54（守 56 上限）。
set -uo pipefail
# 关键：一次性 activate 环境，让 humann 及其子进程(metaphlan/bowtie2/diamond)共享同一 PATH。
# 用逐条 `humann` 会导致 humann 内部调 metaphlan 时子进程 PATH 丢失
# → "CRITICAL ERROR: Can not call software version for metaphlan"（已踩坑，见 plan change-log）。
source /Work_bio/gao/miniforge3/etc/profile.d/conda.sh
conda activate mag_biobakery
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
# metaphlan version-check shim 放到 PATH 最前，修 HUMAnN 3.9 ↔ MetaPhlAn 4.2.4 版本解析崩溃
export PATH="$PROJ/scripts/mpa_shim:$PATH"
# HUMAnN 3.9 只接受 v3/vJun23 的 MetaPhlAn profile；本机装的 metaphlan 默认出 vJan25(被拒)。
# 用 vJun23 db 让内部 metaphlan 产 vJun23 profile → HUMAnN 接受。
MPA_VJ23_DB=/Work_bio/references/Metagenomics/metaphlan_vJun23
MPA_VJ23_IDX=mpa_vJun23_CHOCOPhlAnSGB_202307
export METAPHLAN_DB_DIR=$MPA_VJ23_DB
WORK=$PROJ/work_taxprofiler
CLEAN=$PROJ/clean_reads
HOUT=$PROJ/humann_out
MPA_DB=/Work_bio/references/Metagenomics/metaphlan
MPA_IDX=mpa_vJan25_CHOCOPhlAnSGB_202503
THREADS=18
mkdir -p "$CLEAN" "$HOUT" "$PROJ/logs/humann"

SAMPLES="HFD_AL_4_02_25 HFD_AL_4_03_11 HFD_AL_6_05_12 HFD_AL_6_05_22 HFD_AL_7_06_12 \
         HFD_IF_4_02_25 HFD_IF_4_03_11 HFD_IF_6_05_12 HFD_IF_6_05_22 HFD_IF_7_06_12"

# ---- 1) 准备合并的去宿主 reads（每样本一份，R1+R2 拼一起）----
prep_one(){
  local s=$1
  local out="$CLEAN/${s}.fastq.gz"
  [ -s "$out" ] && { echo "[prep] $s exists"; return 0; }
  local r1 r2
  r1=$(find "$WORK" -name "${s}_${s}_L4.unmapped_1.fastq.gz" 2>/dev/null | head -1)
  r2=$(find "$WORK" -name "${s}_${s}_L4.unmapped_2.fastq.gz" 2>/dev/null | head -1)
  r1=$(readlink -f "$r1"); r2=$(readlink -f "$r2")
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then echo "[prep] MISSING reads for $s"; return 1; fi
  cat "$r1" "$r2" > "$out"
  echo "[prep] $s -> $(du -h $out | cut -f1)"
}
export -f prep_one; export WORK CLEAN
echo "[$(date)] ===== preparing merged host-removed reads ====="
for s in $SAMPLES; do prep_one "$s"; done

# ---- 2) HUMAnN per-sample（3 并行）----
run_one(){
  local s=$1
  local log="$PROJ/logs/humann/${s}.log"
  [ -s "$HOUT/${s}/${s}_pathabundance.tsv" ] && { echo "[humann] $s done"; return 0; }
  echo "[$(date)] [humann] START $s"
  # 内部 metaphlan 用 vJun23 index 产 HUMAnN 兼容的 profile；shim 修版本检查解析
  humann \
     --input "$CLEAN/${s}.fastq.gz" --output "$HOUT/${s}" --output-basename "$s" \
     --threads $THREADS --remove-temp-output \
     --metaphlan-options "-t rel_ab --index $MPA_VJ23_IDX --db_dir $MPA_VJ23_DB --offline --nproc $THREADS" \
     > "$log" 2>&1 && echo "[$(date)] [humann] DONE $s" || echo "[$(date)] [humann] FAIL $s (see $log)"
}
export -f run_one; export PROJ HOUT MPA_VJ23_DB MPA_VJ23_IDX THREADS
echo "[$(date)] ===== running HUMAnN (3 parallel) ====="
printf "%s\n" $SAMPLES | xargs -P 3 -I{} bash -c 'run_one "$@"' _ {}

# ---- 3) 合并 + 标准化 + regroup ----
echo "[$(date)] ===== joining + normalizing tables ====="
MERGE=$PROJ/humann_merged; mkdir -p "$MERGE"
CR=""
for kind in genefamilies pathabundance pathcoverage; do
  humann_join_tables --input "$HOUT" --output "$MERGE/all_${kind}.tsv" --file_name "$kind" --search-subdirectories 2>/dev/null
done
# 相对丰度标准化（CPM 用于 genefamilies，relab 用于 pathways）
humann_renorm_table --input "$MERGE/all_genefamilies.tsv"  --output "$MERGE/all_genefamilies_cpm.tsv"  --units cpm   --update-snames 2>/dev/null
humann_renorm_table --input "$MERGE/all_pathabundance.tsv" --output "$MERGE/all_pathabundance_relab.tsv" --units relab --update-snames 2>/dev/null
# regroup 到 KO / EC（功能通路层面）
humann_regroup_table --input "$MERGE/all_genefamilies_cpm.tsv" --groups uniref90_ko --output "$MERGE/all_ko_cpm.tsv" 2>/dev/null
humann_regroup_table --input "$MERGE/all_genefamilies_cpm.tsv" --groups uniref90_level4ec --output "$MERGE/all_ec_cpm.tsv" 2>/dev/null
echo "[$(date)] ===== HUMANN_ALL_DONE ====="
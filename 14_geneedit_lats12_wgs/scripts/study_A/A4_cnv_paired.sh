#!/bin/bash
# ============================================================================
# Study A / Step 4 — 拷贝数/倍性(配对 tumor vs RO_origin)—— 正式版 2026-07-12
#   两阶段:
#   Part 1+2【现在跑,内存轻~1.5G/样,填闲置CPU】mosdepth 分箱深度 → 相对 origin 的
#            log2 深度比值 → 染色体级增删/非整倍体谱。直接读 CRAM(-f 参考),不转 BAM。
#   Part 3 【精算,随后跑】生成 Control-FREEC 配对 config(control=origin),供审阅后运行。
#   背景:5 个 tumor/编辑细胞为 Brca1/Brca2/Pten CRISPR KO;HRD(Brca1/2缺失)→ 预期
#         拷贝数不稳定/非整倍体。无纯正常对照,以同源 RO_origin 作配对 control。
#   前置:A2 sarek 完成(CRAM 在 output_A/preprocessing/markduplicates/,比对参考=纯GRCm39)。
#   用法: bash A4_cnv_paired.sh              # 跑 Part1+2(mosdepth比值);Part3只生成config
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
REF="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
OUT="$PROJ/analysis_A/cnv"; mkdir -p "$OUT/mosdepth" "$OUT/ratio" "$OUT/freec_configs"
NORMAL="RO_origin"
TUMORS="RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3"
WIN=50000; MAPQ=20; JOBS=3         # JOBS=同时跑的mosdepth数(与B并存,限并发防撞内存)
RUN(){ conda run -n regular_bioinfo "$@"; }
cram(){ ls "$PROJ"/output_A/preprocessing/markduplicates/"$1"/*.cram 2>/dev/null | head -1; }

# ---------- Part 1: mosdepth 分箱深度(限 JOBS 并发) ----------
echo ">> Part 1: mosdepth 分箱深度(${WIN}bp 窗, MAPQ>=$MAPQ, 每样 -t4, 并发 $JOBS)"
run_one(){ local s="$1"; local c; c=$(cram "$s")
  [ -z "$c" ] && { echo "  ERROR: 无 $s CRAM"; return 1; }
  RUN mosdepth --by "$WIN" -Q "$MAPQ" -n --fast-mode -t 4 -f "$REF" "$OUT/mosdepth/$s" "$c" \
    && echo "  ✔ mosdepth $s"; }
i=0
for s in $NORMAL $TUMORS; do
  run_one "$s" &
  i=$((i+1)); [ $((i % JOBS)) -eq 0 ] && wait
done
wait
echo "  Part1 完成。"

# ---------- Part 2: log2 深度比值 + 染色体级汇总 ----------
echo ">> Part 2: 计算各 tumor 相对 $NORMAL 的 log2 深度比值(常染色体, 自身中位归一)"
RUN python3 - "$OUT" "$NORMAL" "$TUMORS" <<'PY'
import sys, gzip, math, statistics as st
from collections import defaultdict
OUT, NORMAL, TUMORS = sys.argv[1], sys.argv[2], sys.argv[3].split()
AUTO = {f"chr{i}" for i in range(1,20)}
def load(s):
    w={}
    with gzip.open(f"{OUT}/mosdepth/{s}.regions.bed.gz","rt") as f:
        for l in f:
            c,a,b,d=l.split()[:4]
            if c in AUTO: w[(c,int(a))]=float(d)
    return w
def med_nonzero(w):
    v=[x for x in w.values() if x>0]; return st.median(v) if v else 0
norm=load(NORMAL); nmed=med_nonzero(norm)
print(f"  {NORMAL} 常染色体中位深度={nmed:.1f}x  窗口数={len(norm)}")
for t in TUMORS:
    tw=load(t); tmed=med_nonzero(tw)
    rows=[]; devwin=0; tot=0; chrom_l2=defaultdict(list)
    for k in sorted(set(tw)&set(norm)):
        nd,td=norm[k],tw[k]
        if nd<=0 or td<=0 or nmed<=0 or tmed<=0: continue
        l2=math.log2((td/tmed)/(nd/nmed))
        rows.append((k[0],k[1],k[1]+50000,round(l2,3)))
        chrom_l2[k[0]].append(l2); tot+=1
        if abs(l2)>0.3: devwin+=1
    with open(f"{OUT}/ratio/{t}.log2ratio.bed","w") as o:
        o.write("chrom\tstart\tend\tlog2ratio_vs_origin\n")
        for r in rows: o.write("\t".join(map(str,r))+"\n")
    aneu = 100.0*devwin/tot if tot else 0
    gains=[c for c in chrom_l2 if st.median(chrom_l2[c])>0.3]
    losses=[c for c in chrom_l2 if st.median(chrom_l2[c])<-0.3]
    with open(f"{OUT}/ratio/{t}.chrom_summary.tsv","w") as o:
        o.write("chrom\tmedian_log2ratio\tcall\n")
        for c in sorted(chrom_l2, key=lambda x:int(x[3:])):
            m=st.median(chrom_l2[c]); call="gain" if m>0.3 else "loss" if m<-0.3 else "neutral"
            o.write(f"{c}\t{m:.3f}\t{call}\n")
    print(f"  {t}: 中位深度={tmed:.1f}x 非整倍体分数(|log2|>0.3窗占比)={aneu:.1f}%  "
          f"臂级增={gains or '无'} 减={losses or '无'}")
PY
echo "  Part2 完成 → $OUT/ratio/"

# ---------- Part 3: 生成 Control-FREEC 配对 config(供审阅,精算随后跑) ----------
echo ">> Part 3: 生成 FREEC 配对 config(不运行,待审阅)"
FAI="$REF.fai"; NCRAM=$(cram "$NORMAL")
for t in $TUMORS; do
  TCRAM=$(cram "$t")
  cat > "$OUT/freec_configs/${t}_vs_${NORMAL}.freec.config" <<CFG
[general]
chrLenFile = $FAI
ploidy = 2
window = $WIN
outputDir = $OUT/freec/${t}
maxThreads = 6
samtools = $(RUN which samtools 2>/dev/null || echo samtools)
[sample]
mateFile = $TCRAM
inputFormat = CRAM
mateOrientation = FR
[control]
mateFile = $NCRAM
inputFormat = CRAM
mateOrientation = FR
# ⚠ CRAM 输入需 FREEC 支持 htslib;若不支持则先转 BAM(用 output_A CRAM + $REF)。审阅时确认。
CFG
done
echo "  已生成 5 个 FREEC config → $OUT/freec_configs/(审阅后再跑 FREEC 精算)"
echo "DONE A4: mosdepth 比值谱已出(analysis_A/cnv/ratio/);FREEC config 待审阅运行。"

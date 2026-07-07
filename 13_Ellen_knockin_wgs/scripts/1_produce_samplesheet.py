#!/usr/bin/env python3
# ============================================================================
# Step 1 — 生成 nf-core/sarek samplesheet
#   扫描 Ellen fastq 目录，从文件名解析样本名(下划线前缀=打靶品系)，
#   输出：
#     samplesheet_full.csv         全部 6 样
#     samplesheet_trial_RAGH.csv   仅 RAGH_153（管线首跑验证用）
#   sarek 列：patient,sample,lane,fastq_1,fastq_2
#   （每样各自为一个 patient；纯 germline；无 tumor-normal 配对）
# ============================================================================
import os, re, glob, csv, sys

FASTQ_DIR = "/home/gao/Dropbox/Ellen"
OUT_DIR   = os.path.dirname(os.path.abspath(__file__))

# 文件名: 23KCYFLT4_7_0469165304_<SAMPLE>_<Snn>_L007_R1_001.fastq.gz
#   <SAMPLE> 如 RAGH_153 / MTTH_284 / CD1A_B125
PAT = re.compile(r"_(?P<name>[A-Za-z0-9]+_[A-Za-z0-9]+)_(?P<snum>S\d+)_L(?P<lane>\d+)_R1_001\.fastq\.gz$")

rows = []
for r1 in sorted(glob.glob(os.path.join(FASTQ_DIR, "*_R1_001.fastq.gz"))):
    m = PAT.search(os.path.basename(r1))
    if not m:
        print("WARN: 无法解析文件名，跳过:", r1); continue
    name = m.group("name")
    lane = "L%03d" % int(m.group("lane"))
    r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    if not os.path.exists(r2):
        print("ERROR: 缺少 R2:", r2); sys.exit(1)
    rows.append({"patient": name, "sample": name, "lane": lane,
                 "fastq_1": r1, "fastq_2": r2})

if not rows:
    print("ERROR: 未找到任何 fastq"); sys.exit(1)

FIELDS = ["patient", "sample", "lane", "fastq_1", "fastq_2"]

def write(path, subset):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader()
        for r in subset: w.writerow(r)
    print(f"  写入 {path}  ({len(subset)} 样)")

write(os.path.join(OUT_DIR, "samplesheet_full.csv"), rows)
trial = [r for r in rows if r["sample"] == "RAGH_153"]
if not trial:
    print("ERROR: 试跑目标 RAGH_153 未找到"); sys.exit(1)
write(os.path.join(OUT_DIR, "samplesheet_trial_RAGH.csv"), trial)

print("\n确认表（全部样本）：")
print(f"  {'sample':<12} {'lane':<6} R1")
for r in rows:
    print(f"  {r['sample']:<12} {r['lane']:<6} {os.path.basename(r['fastq_1'])}")

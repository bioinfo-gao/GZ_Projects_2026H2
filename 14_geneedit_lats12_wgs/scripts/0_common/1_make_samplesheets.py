#!/usr/bin/env python3
# ============================================================================
# Step 1 — 由 ../docs/sample_info.tsv 生成 sarek samplesheet（两份）
#   A_somatic.csv  : patient=RO, RO_origin=normal(status0), 其余=tumor(status1)
#   B_germline.csv : 每样各自 patient, status0（germline）
#   sarek 列: patient,sample,status,lane,fastq_1,fastq_2
#   fastq: /home/gao/Dropbox/JinPeng/*_<fastq_prefix>_L006_R{1,2}_001.fastq.gz
# ============================================================================
import csv, glob, os, sys

TSV = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs/docs/sample_info.tsv"
FASTQ_DIR = "/home/gao/Dropbox/JinPeng"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
LANE = "L006"

def find_fq(prefix, r):
    hits = glob.glob(os.path.join(FASTQ_DIR, f"*_{prefix}_{LANE}_R{r}_001.fastq.gz"))
    if len(hits) != 1:
        sys.exit(f"ERROR: {prefix} R{r} 找到 {len(hits)} 个文件（应为1）")
    return hits[0]

rows = list(csv.DictReader(open(TSV), delimiter="\t"))
A, B = [], []
for x in rows:
    fq1, fq2 = find_fq(x["fastq_prefix"], 1), find_fq(x["fastq_prefix"], 2)
    if x["study"] == "A":
        A.append(dict(patient="RO", sample=x["sample_id"], status=x["status"],
                      lane=LANE, fastq_1=fq1, fastq_2=fq2))
    else:  # B germline: 每样自成 patient
        B.append(dict(patient=x["sample_id"], sample=x["sample_id"], status="0",
                      lane=LANE, fastq_1=fq1, fastq_2=fq2))

FIELDS = ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
def write(path, data):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(data)
    print(f"  写入 {path}（{len(data)} 样）")

write(os.path.join(OUT_DIR, "A_somatic.csv"), A)
write(os.path.join(OUT_DIR, "B_germline.csv"), B)

print("\nStudy A（patient=RO；status 0=normal/origin, 1=tumor）:")
for r in A: print(f"  {r['sample']:<12} status={r['status']}  {os.path.basename(r['fastq_1'])}")
print("Study B（germline）:")
for r in B: print(f"  {r['sample']:<12} {os.path.basename(r['fastq_1'])}")

#!/usr/bin/env python3
"""
生成本项目所有 pipeline 输入表：
  - taxprofiler:  samplesheet_taxprofiler.csv  (sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta)
  - taxprofiler:  databases.csv                (kraken2 + bracken + metaphlan，指向本地共享库)
  - nf-core/mag:  samplesheet_mag.csv          (sample,group,short_reads_1,short_reads_2,long_reads / 5.4.2 schema)
数据源：/home/gao/Dropbox/QTE_26_06_25_001_Daniel_Mendes/HFD_*  (PE150 短读长，mouse stool)
"""
import glob, os, csv

DATA = "/home/gao/Dropbox/QTE_26_06_25_001_Daniel_Mendes"
OUT  = os.path.dirname(os.path.abspath(__file__)) + "/.."
PLATFORM = "ILLUMINA"

samples = []
for d in sorted(glob.glob(f"{DATA}/HFD_*")):
    name = os.path.basename(d)
    r1 = sorted(glob.glob(f"{d}/*_1.fq.gz"))
    r2 = sorted(glob.glob(f"{d}/*_2.fq.gz"))
    assert len(r1) == 1 and len(r2) == 1, f"{name}: expected 1 R1/R2 pair, got {r1} {r2}"
    group = "AL" if "_AL_" in name else "IF"
    samples.append(dict(sample=name, group=group, r1=r1[0], r2=r2[0]))

assert len(samples) == 10, f"expected 10 samples, got {len(samples)}"

# ---- taxprofiler samplesheet ----
with open(f"{OUT}/samplesheet_taxprofiler.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "run_accession", "instrument_platform", "fastq_1", "fastq_2", "fasta"])
    for s in samples:
        w.writerow([s["sample"], f"{s['sample']}_L4", PLATFORM, s["r1"], s["r2"], ""])

# ---- taxprofiler databases.csv (读长 PE150 -> bracken -r 150) ----
with open(f"{OUT}/databases.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tool", "db_name", "db_params", "db_type", "db_path"])
    k2 = "/Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226/k2_standard_08_GB_20260226.tar.gz"
    mpa = "/Work_bio/references/Metagenomics/metaphlan/"
    # db_name MUST be unique per row — 共用同一 db_name 会让 Bracken 跑两遍并在
    # BRACKEN_COMBINEBRACKENOUTPUTS 触发 input file name collision（已踩坑，2026-07-18）。
    w.writerow(["kraken2",   "k2s8_kraken2", "",        "short", k2])
    w.writerow(["bracken",   "k2s8_bracken", ";-r 150", "short", k2])
    w.writerow(["metaphlan", "mpa_vJan25",   "",        "short", mpa])

# ---- nf-core/mag samplesheet (5.4.2) ----
with open(f"{OUT}/samplesheet_mag.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "group", "short_reads_1", "short_reads_2", "long_reads"])
    for s in samples:
        w.writerow([s["sample"], s["group"], s["r1"], s["r2"], ""])

print(f"OK: {len(samples)} samples")
for s in samples:
    print(f"  {s['sample']:<18} group={s['group']}")

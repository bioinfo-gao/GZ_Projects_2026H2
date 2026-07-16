#!/usr/bin/env python3
"""Generate the nf-core/sarek samplesheet.csv from the canonical sample_info.tsv.

sarek schema (fixed): patient,sample,sex,status,lane,fastq_1,fastq_2
  - germline: each sample is its own `patient`, status=0 (no tumor/normal pairing)
  - sex=NA  -> sarek estimates from coverage; we don't assume
  - lane    -> read-group identifier (single lane L4 here)
Client-facing names live only in sample_info.tsv (single source of truth).
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
info = ROOT / "sample_info.tsv"
out = ROOT / "scripts" / "samplesheet.csv"

rows = list(csv.DictReader(info.open(), delimiter="\t"))
with out.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["patient", "sample", "sex", "status", "lane", "fastq_1", "fastq_2"])
    for r in rows:
        w.writerow([r["sample"], r["sample"], "NA", "0",
                    r["flowcell_lane"], r["fastq_1"], r["fastq_2"]])
print(f"wrote {out} ({len(rows)} samples)")

#!/usr/bin/env python3
"""
P16 Wenliang addendum -- gene-overlap annotation for SV (Manta, PASS) and CNV (CNVkit,
non-diploid) calls. Lightweight bedtools+GTF overlap, NOT a full AnnotSV clinical
annotation -- matches the ~1h scope quoted to the client for this item.
env: regular_bioinfo (pysam, bedtools on PATH)
"""
import csv
import re
import subprocess
from pathlib import Path

import pysam

PROJ = Path("/home/gao/projects_2026H2/16_Wenliang_Pan_wgs")
SRC = PROJ / "custom_research_report_20260717" / "structural_cnv"
OUT = PROJ / "custom_research_report_20260720_addendum" / "structural_cnv_annotation"
GTF = Path("/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf")
CACHE_DIR = PROJ / "scripts" / "rds_cache_addendum"
GENEBED = CACHE_DIR / "genes.grch38.bed"

OUT.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ---- 1. gene BED (chrom, start0, end, gene_name), sorted ----
if not GENEBED.exists():
    gene_re = re.compile(r'gene_name "([^"]+)"')
    rows = []
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            f_ = line.rstrip("\n").split("\t")
            if f_[2] != "gene":
                continue
            m = gene_re.search(f_[8])
            gene = m.group(1) if m else "NA"
            rows.append((f_[0], int(f_[3]) - 1, int(f_[4]), gene))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(GENEBED, "w") as f:
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
print(f"gene BED: {sum(1 for _ in open(GENEBED))} genes")


def collapse_overlap(bed_a: Path, bed_b: Path, extra_cols: int):
    """bedtools intersect -a A -b B -wa -wb, collapse B's gene col (last field) per A record."""
    out = subprocess.run(
        ["bedtools", "intersect", "-a", str(bed_a), "-b", str(bed_b), "-wa", "-wb"],
        check=True, capture_output=True, text=True,
    ).stdout
    genes = {}
    for line in out.splitlines():
        f_ = line.split("\t")
        key = tuple(f_[:3 + extra_cols])
        gene = f_[-1]
        genes.setdefault(key, [])
        if gene not in genes[key]:
            genes[key].append(gene)
    return genes


for sample in ("Sample_A", "Sample_B"):
    print(f"=== {sample} ===")

    # ---- Manta SV (PASS only) -> BED ----
    manta_vcf = SRC / f"{sample}.manta.diploid_sv.vcf.gz"
    sv_bed = OUT / f"{sample}.manta_pass.bed"
    sv_records = []
    with pysam.VariantFile(str(manta_vcf)) as vf:
        for rec in vf:
            if list(rec.filter) != ["PASS"]:
                continue
            end = rec.info.get("END", rec.pos + 1)
            if end <= rec.pos:
                end = rec.pos + 1
            svtype = rec.info.get("SVTYPE", "NA")
            svlen = rec.info.get("SVLEN", "NA")
            if isinstance(svlen, tuple):
                svlen = svlen[0]
            sv_records.append((rec.chrom, rec.pos - 1, end, rec.id, svtype, str(svlen)))
    sv_records.sort(key=lambda r: (r[0], r[1]))
    with open(sv_bed, "w") as f:
        for r in sv_records:
            f.write("\t".join(str(x) for x in r) + "\n")
    n_sv = len(sv_records)
    print(f"  Manta PASS SVs: {n_sv}")

    gene_hits = collapse_overlap(sv_bed, GENEBED, extra_cols=3)
    out_tsv = OUT / f"{sample}.manta_pass.gene_overlap.tsv"
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "start", "end", "sv_id", "sv_type", "sv_len", "genes_overlapping"])
        for r in sv_records:
            key = tuple(str(x) for x in r[:3]) + (r[3], r[4], r[5])
            key_lookup = tuple(str(x) for x in r[:3])
            # collapse_overlap keyed on (chrom,start,end,id,type,len) since extra_cols=3
            full_key = (str(r[0]), str(r[1]), str(r[2]), r[3], r[4], r[5])
            genes = gene_hits.get(full_key)
            if genes:
                w.writerow([r[0], r[1], r[2], r[3], r[4], r[5], ",".join(genes)])
    n_hit = sum(1 for _ in open(out_tsv)) - 1
    print(f"  Manta SVs overlapping >=1 gene: {n_hit} / {n_sv}")
    sv_bed.unlink()

    # ---- CNVkit non-diploid segments (cn != 2) -> BED ----
    cns = SRC / f"{sample}.md.call.cns"
    cnv_bed = OUT / f"{sample}.cnv_nondiploid.bed"
    cnv_records = []
    with open(cns) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            f_ = line.rstrip("\n").split("\t")
            cn = int(f_[idx["cn"]])
            if cn == 2:
                continue
            cnv_records.append((f_[idx["chromosome"]], int(f_[idx["start"]]), int(f_[idx["end"]]),
                                 f_[idx["log2"]], str(cn)))
    cnv_records.sort(key=lambda r: (r[0], r[1]))
    with open(cnv_bed, "w") as f:
        for r in cnv_records:
            f.write("\t".join(str(x) for x in r) + "\n")
    n_cnv = len(cnv_records)
    print(f"  CNVkit non-diploid segments (cn!=2): {n_cnv}")

    gene_hits = collapse_overlap(cnv_bed, GENEBED, extra_cols=2)
    out_tsv = OUT / f"{sample}.cnv_nondiploid.gene_overlap.tsv"
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "start", "end", "log2", "cn", "genes_overlapping"])
        for r in cnv_records:
            full_key = (str(r[0]), str(r[1]), str(r[2]), r[3], r[4])
            genes = gene_hits.get(full_key)
            if genes:
                w.writerow([r[0], r[1], r[2], r[3], r[4], ",".join(genes)])
    n_cnvhit = sum(1 for _ in open(out_tsv)) - 1
    print(f"  CNV segments overlapping >=1 gene: {n_cnvhit} / {n_cnv}")
    cnv_bed.unlink()

print("DONE 11_sv_cnv_gene_overlap")

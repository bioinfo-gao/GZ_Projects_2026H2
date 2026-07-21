#!/usr/bin/env python3
"""
P16 Wenliang addendum -- two client-requested additions, both derived from the existing
prioritised (rare + functional/ClinVar) variant lists, NOT the raw ~4-5M genome-wide callset:

  1. Zygosity (het/hom) column added to a NEW copy of all_samples.prioritised.tsv
     (the original 0717 delivery is left untouched -- see explanation doc for why).
  2. Shared vs. sample-specific variant comparison between Sample_A and Sample_B,
     restricted to the prioritised list (comparing the full raw callset between two
     unrelated-sex individuals would mostly surface common population polymorphisms,
     not anything actionable).

Requires: pysam. Run: conda run -n regular_bioinfo python scripts/12_zygosity_and_comparison.py
"""
import csv
import gzip
from pathlib import Path

PROJ = Path("/home/gao/projects_2026H2/16_Wenliang_Pan_wgs")
PRIORITISED = PROJ / "custom_research_report_20260717" / "annotation_prioritised" / "all_samples.prioritised.tsv"
VCF_DIR = PROJ / "output_results" / "annotation_gnomad_clinvar"
OUT_ZYG = PROJ / "custom_research_report_20260720_addendum" / "annotation_prioritised_updated" / "all_samples.prioritised_with_zygosity.tsv"
OUT_CMP_DIR = PROJ / "custom_research_report_20260720_addendum" / "comparison_A_vs_B"
OUT_CMP_DIR.mkdir(parents=True, exist_ok=True)
OUT_ZYG.parent.mkdir(parents=True, exist_ok=True)

def zyg_from_gt(gt):
    """Generic GT -> zygosity, handling multi-allelic GTs (e.g. 1/2) correctly."""
    if not gt or gt in (".", "./.", ".|."):
        return "no_call"
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    alleles = gt.split(sep) if sep else [gt]
    if len(set(alleles)) == 1:
        return "homozygous_ref" if alleles[0] == "0" else "homozygous_alt"
    return "heterozygous"


def load_gt_lookup(sample):
    """(chrom,pos,ref,alt) -> GT, read from the sample's annotated single-sample VCF.
    ALT is kept EXACTLY as it appears in the VCF (comma-joined for multi-allelic sites),
    matching how the upstream prioritised.tsv (script 5) stores it -- do not split on comma."""
    vcf = VCF_DIR / f"{sample}.gnomad.clinvar.vcf.gz"
    lut = {}
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = c[0], c[1], c[2], c[3], c[4]
            fmt = c[8].split(":")
            sample_fields = c[9].split(":")
            gt = sample_fields[fmt.index("GT")] if "GT" in fmt else ""
            lut[(chrom, pos, ref, alt)] = gt
    return lut


# ---- 1. load prioritised table, add zygosity ----
with open(PRIORITISED) as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames

gt_lut = {s: load_gt_lookup(s) for s in ("Sample_A", "Sample_B")}

for r in rows:
    gt = gt_lut[r["sample"]].get((r["chrom"], r["pos"], r["ref"], r["alt"]), "")
    r["zygosity"] = zyg_from_gt(gt)

new_fieldnames = fieldnames + ["zygosity"]
with open(OUT_ZYG, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=new_fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"zygosity table: {len(rows)} rows -> {OUT_ZYG}")
n_missing = sum(1 for r in rows if r["zygosity"] in ("", "NA"))
print(f"  rows with no GT match (unexpected if >0): {n_missing}")

# ---- 2. shared vs sample-specific comparison, on the prioritised list ----
by_sample = {"Sample_A": {}, "Sample_B": {}}
for r in rows:
    key = (r["chrom"], r["pos"], r["ref"], r["alt"])
    by_sample[r["sample"]][key] = r

keys_a = set(by_sample["Sample_A"])
keys_b = set(by_sample["Sample_B"])
shared = keys_a & keys_b
only_a = keys_a - keys_b
only_b = keys_b - keys_a

cols = ["chrom", "pos", "ref", "alt", "gene", "consequence", "clinvar_sig", "clinvar_disease",
        "status", "zygosity_Sample_A", "zygosity_Sample_B"]


def row_for(key, status):
    ra = by_sample["Sample_A"].get(key)
    rb = by_sample["Sample_B"].get(key)
    base = ra or rb
    return {
        "chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
        "gene": base["gene"], "consequence": base["consequence"],
        "clinvar_sig": base["clinvar_sig"], "clinvar_disease": base["clinvar_disease"],
        "status": status,
        "zygosity_Sample_A": ra["zygosity"] if ra else "absent",
        "zygosity_Sample_B": rb["zygosity"] if rb else "absent",
    }


out_rows = (
    [row_for(k, "shared") for k in sorted(shared)]
    + [row_for(k, "unique_to_Sample_A") for k in sorted(only_a)]
    + [row_for(k, "unique_to_Sample_B") for k in sorted(only_b)]
)
out_tsv = OUT_CMP_DIR / "prioritised_variant_comparison_A_vs_B.tsv"
with open(out_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)

print(f"comparison table: {len(out_rows)} rows -> {out_tsv}")
print(f"  shared: {len(shared)}  unique_to_Sample_A: {len(only_a)}  unique_to_Sample_B: {len(only_b)}")
print(f"  (prioritised list totals: Sample_A={len(keys_a)}, Sample_B={len(keys_b)})")

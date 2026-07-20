#!/usr/bin/env python3
# P17 MAG — build the non-redundant final MAG catalog.
# nf-core/mag QCs (CheckM2/BUSCO/QUAST) and GTDB-Tk-classifies the raw per-binner bins
# (MetaBAT2/MaxBin2/SemiBin2, 1060 rows incl. cross-binner duplicates of the same genome),
# while DASTool separately selects a non-redundant 184-bin consensus set per group (renamed
# "<Binner>Refined-group-XX.NNN[_sub].fa"). This script maps each DASTool-refined bin back to
# its originating raw bin (by binner+group+number, stripping "Refined"/"_sub") to attach that
# bin's QC + taxonomy, giving one row per genuinely distinct recovered genome.
import csv, os, re

PROJ = "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
MAG  = os.path.join(PROJ, "output_results_mag")
OUT  = os.path.join(PROJ, "custom_research_report_20260720", "assembly_binning")

with open(os.path.join(MAG, "GenomeBinning", "bin_summary.tsv")) as f:
    rows = {r["bin"]: r for r in csv.DictReader(f, delimiter="\t")}

def map_name(n):
    base = n.replace("Refined", "")
    base = re.sub(r"_sub(?=\.fa$)", "", base)
    if "SemiBin2" in base:
        base = re.sub(r"(group-(?:AL|IF))\.(\d+)\.fa$", r"\1_\2.fa", base)
    return base

def quality_tier(comp, cont):
    if comp is None or cont is None:
        return "NA"
    if comp >= 90 and cont < 5:
        return "High"
    if comp >= 50 and cont < 10:
        return "Medium"
    return "Low"

dastool_dir = os.path.join(MAG, "GenomeBinning", "DASTool", "bins")
out_rows, unmapped = [], []
for fname in sorted(os.listdir(dastool_dir)):
    group = "AL" if "group-AL" in fname else "IF"
    binner = re.match(r"MEGAHIT-(\w+?)Refined-", fname).group(1)
    raw_name = map_name(fname)
    src = rows.get(raw_name)
    if src is None:
        unmapped.append(fname)
        continue
    comp = float(src["Completeness_checkm2"]) if src["Completeness_checkm2"] not in ("", "NA") else None
    cont = float(src["Contamination_checkm2"]) if src["Contamination_checkm2"] not in ("", "NA") else None
    out_rows.append({
        "final_bin": fname,
        "group": group,
        "origin_binner": binner,
        "raw_bin_matched": raw_name,
        "completeness_pct": comp,
        "contamination_pct": cont,
        "quality_tier_MIMAG": quality_tier(comp, cont),
        "genome_size_bp": src["Genome_Size_checkm2"],
        "n50_bp": src["N50_quast"],
        "gc_pct": src["GC (%)_quast"],
        "num_contigs": src["# contigs_quast"],
        "gtdb_classification": src["classification_gtdbtk"],
    })

assert not unmapped, f"unmapped DASTool bins: {unmapped}"
assert len(out_rows) == 184, f"expected 184 final bins, got {len(out_rows)}"

os.makedirs(OUT, exist_ok=True)
out_path = os.path.join(OUT, "final_MAG_catalog.tsv")
fields = list(out_rows[0].keys())
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)

# summary table: group x quality tier counts + classified fraction
from collections import Counter
summary = []
for g in ("AL", "IF"):
    grows = [r for r in out_rows if r["group"] == g]
    tiers = Counter(r["quality_tier_MIMAG"] for r in grows)
    classified = sum(1 for r in grows if r["gtdb_classification"].strip())
    summary.append({
        "group": g,
        "total_MAGs": len(grows),
        "high_quality_gt90comp_lt5cont": tiers.get("High", 0),
        "medium_quality_50to90comp_lt10cont": tiers.get("Medium", 0),
        "low_quality": tiers.get("Low", 0),
        "gtdb_classified": classified,
        "gtdb_unclassified": len(grows) - classified,
    })
sum_path = os.path.join(OUT, "mag_quality_summary_by_group.tsv")
with open(sum_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(summary)

print(f"final_MAG_catalog.tsv: {len(out_rows)} rows -> {out_path}")
for s in summary:
    print(s)

#!/usr/bin/env python3
"""Project 16 — identify rare + potentially functional germline variants.

Input : output_results/annotation_gnomad_clinvar/<sample>.gnomad.clinvar.vcf.gz
        (VEP CSQ from sarek + gnomAD_AF + CLNSIG/CLNDN from step 4)
Output: per-sample TSV of prioritised variants + a combined summary.

Prioritisation logic (report in plan §4):
  RARE       : gnomAD_AF < 0.001 OR absent from gnomAD
  FUNCTIONAL : VEP consequence is HIGH/MODERATE impact (LoF, missense, splice, ...)
               OR ClinVar CLNSIG contains Pathogenic / Likely_pathogenic
  PASS filter: FILTER == PASS (or '.')
A variant is flagged if RARE and (FUNCTIONAL). ClinVar P/LP is always surfaced
regardless of frequency (known clinically-actionable).

Requires: pysam (regular_bioinfo). Run: conda run -n regular_bioinfo python scripts/5_rare_functional_filter.py
"""
import glob
import gzip
import os
import re

PROJ = "/home/gao/projects_2026H2/16_Wenliang_Pan_wgs"
INDIR = f"{PROJ}/output_results/annotation_gnomad_clinvar"
OUTDIR = f"{PROJ}/output_results/prioritised_variants"
os.makedirs(OUTDIR, exist_ok=True)

RARE_AF = 1e-3
HIGH_MOD = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
}
PATHO = re.compile(r"[Pp]athogenic")


def parse_csq_format(header_lines):
    for h in header_lines:
        if h.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format: ([^\">]+)", h)
            if m:
                return m.group(1).split("|")
    return []


def info_dict(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def process(vcf):
    sample = os.path.basename(vcf).split(".")[0]
    header, rows = [], []
    op = gzip.open(vcf, "rt")
    csq_fields = []
    n_total = n_flag = 0
    for line in op:
        if line.startswith("##"):
            header.append(line)
            continue
        if line.startswith("#"):
            csq_fields = parse_csq_format(header)
            continue
        n_total += 1
        c = line.rstrip("\n").split("\t")
        chrom, pos, _id, ref, alt, _q, filt, info = c[:8]
        if filt not in ("PASS", "."):
            continue
        d = info_dict(info)
        try:
            af = float(d.get("gnomAD_AF", "nan"))
        except ValueError:
            af = float("nan")
        rare = (af != af) or (af < RARE_AF)  # nan -> absent -> rare
        clnsig = d.get("CLNSIG", "")
        clinvar_plp = bool(PATHO.search(clnsig))
        # worst VEP consequence + gene across transcripts
        cons, genes = set(), set()
        if "CSQ" in d and csq_fields:
            gi = csq_fields.index("Consequence") if "Consequence" in csq_fields else None
            sy = csq_fields.index("SYMBOL") if "SYMBOL" in csq_fields else None
            for tx in d["CSQ"].split(","):
                parts = tx.split("|")
                if gi is not None and gi < len(parts):
                    cons.update(parts[gi].split("&"))
                if sy is not None and sy < len(parts) and parts[sy]:
                    genes.add(parts[sy])
        functional = bool(cons & HIGH_MOD)
        keep = clinvar_plp or (rare and functional)
        if not keep:
            continue
        n_flag += 1
        rows.append([
            sample, chrom, pos, ref, alt,
            ("absent" if af != af else f"{af:.3g}"),
            ";".join(sorted(genes)) or ".",
            ";".join(sorted(cons & HIGH_MOD)) or ";".join(sorted(cons)) or ".",
            clnsig or ".", d.get("CLNDN", ".").replace("_", " "),
            "ClinVar_P/LP" if clinvar_plp else "rare_functional",
        ])
    op.close()
    out = f"{OUTDIR}/{sample}.prioritised.tsv"
    cols = ["sample", "chrom", "pos", "ref", "alt", "gnomAD_AF", "gene",
            "consequence", "clinvar_sig", "clinvar_disease", "flag"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{sample}: {n_total} PASS-eval'd -> {n_flag} prioritised -> {out}")
    return rows, cols


def main():
    vcfs = sorted(glob.glob(f"{INDIR}/*.gnomad.clinvar.vcf.gz"))
    if not vcfs:
        raise SystemExit(f"No annotated VCFs in {INDIR}; run step 4 first.")
    allrows, cols = [], None
    for v in vcfs:
        rows, cols = process(v)
        allrows += rows
    comb = f"{OUTDIR}/all_samples.prioritised.tsv"
    with open(comb, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in allrows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"combined: {len(allrows)} rows -> {comb}")


if __name__ == "__main__":
    main()

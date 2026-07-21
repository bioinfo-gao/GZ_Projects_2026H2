# Addendum to Wenliang_Pan_WGS_0717 — Comparison, Zygosity, and SV/CNV Annotation

**Date:** 2026-07-20
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Species:** *Homo sapiens* (GRCh38) | **Sample/Tissue:** Sample_A, Sample_B (genomic DNA, tissue not specified by client)

> **Scope of this addendum.** This folder answers the three follow-up questions raised on 2026-07-20
> regarding the `custom_research_report_20260717/` delivery. It is a **supplement, not a
> re-analysis** — the original delivery folder and every file in it are unchanged. Nothing here
> alters the conclusions already reported on 2026-07-17; it adds detail that was requested afterward.

## 1. Comparison analysis: shared vs. sample-specific variants (Sample_A vs Sample_B)

**File:** `comparison_A_vs_B/prioritised_variant_comparison_A_vs_B.tsv`

**What was compared.** This comparison is restricted to the **prioritised (rare + functional /
ClinVar) variant list** already delivered in `annotation_prioritised/all_samples.prioritised.tsv`
— it is **not** a comparison of the full raw genome-wide callset (~4-5 million variants/sample).
Comparing the raw callsets directly would be dominated by common population polymorphisms shared
by essentially all humans and would not answer a scientifically useful question; restricting to
the already-curated rare/functional list is what makes "shared" and "unique" mutations meaningful.
If you specifically need the full raw-callset overlap (not restricted to rare/functional variants),
let us know and we can generate that separately.

**Result:**

| | Sample_A | Sample_B |
| :--- | :---: | :---: |
| Prioritised variants (total) | 4,294 | 4,342 |
| **Shared** with the other sample | 2,144 | 2,144 |
| **Unique** to this sample | 2,150 | 2,198 |

**Important context — genetic sex.** As already reported in §6 of the main report, Sample_A is
genetically male (46,XY) and Sample_B is genetically female (46,XX): they are two distinct
individuals, not a matched pair from the same person. The "shared" variants here are therefore
autosomal variants both individuals happen to carry (a mix of coincidence and any true relatedness
between them), not evidence of common origin by itself. If Sample_A and Sample_B are meant to be
biological relatives, this table is still the right comparison — just note that a large "shared"
count is expected for any two humans at the rare-variant level given a large enough prioritised list,
and does not by itself quantify degree of relatedness (which would require a dedicated kinship
analysis — available on request, additional scope).

Each row also carries per-sample **zygosity** (see §2) so you can see, for a shared variant, whether
each person carries it in one copy or two.

## 2. Zygosity (heterozygous / homozygous)

**File:** `annotation_prioritised_updated/all_samples.prioritised_with_zygosity.tsv`

Identical to the original `all_samples.prioritised.tsv` (same 8,636 rows, same columns), with one
column appended: **`zygosity`**, derived from the standard VCF genotype (`GT`) field in the
annotated callset:

- `heterozygous` — one reference + one alternate allele (`0/1`), or two different alternate
  alleles at a multi-allelic site (`1/2`)
- `homozygous_alt` — both alleles are the same alternate allele (`1/1`)

7,565 of 8,636 prioritised variants are heterozygous, 1,071 are homozygous-alt — consistent with the
expected genome-wide het:hom ratio for germline variants.

## 3. Structural variant (SV) and copy-number variant (CNV) gene annotation

**Files:** `structural_cnv_annotation/{Sample_A,Sample_B}.manta_pass.gene_overlap.tsv`,
`structural_cnv_annotation/{Sample_A,Sample_B}.cnv_nondiploid.gene_overlap.tsv`

The original delivery included only the raw Manta/TIDDIT/CNVkit callsets (VCFs and a `.cns`
segment file), with no gene-level annotation — this addendum adds that layer:

- **SV (Manta, PASS calls only):** each SV's genomic interval was intersected against GENCODE v45
  gene models (GRCh38) to list overlapping gene(s). 3,815/6,349 (Sample_A) and 4,367/7,209
  (Sample_B) PASS SVs overlap at least one gene.
- **CNV (CNVkit):** restricted to **non-diploid segments only** (copy number ≠ 2 — the segments
  that are potentially biologically interesting; diploid/CN=2 background segments are not
  annotated). 142/192 (Sample_A) and 123/162 (Sample_B) non-diploid segments overlap at least one
  gene.

**What this is, and isn't.** This is a lightweight gene-overlap annotation (bedtools + GENCODE),
matching the ~1 hour of work quoted for this item — it lists which gene(s) each SV/CNV call
physically overlaps. It is **not** a clinical-grade annotation (no pathogenicity classification,
dosage-sensitivity, or curated disease-gene database cross-reference, the way VEP + gnomAD +
ClinVar work for the SNV/indel calls). If you need that deeper level of SV/CNV interpretation
(e.g., AnnotSV-style pathogenicity flags), that is additional scope beyond this addendum — happy
to quote it separately.

## Deliverable Files

```
custom_research_report_20260720_addendum/
├── README_addendum_0720.md                              <- this document
├── comparison_A_vs_B/
│   └── prioritised_variant_comparison_A_vs_B.tsv         <- shared / unique_to_Sample_A / unique_to_Sample_B
├── annotation_prioritised_updated/
│   └── all_samples.prioritised_with_zygosity.tsv         <- original table + zygosity column
└── structural_cnv_annotation/
    ├── Sample_A.manta_pass.gene_overlap.tsv
    ├── Sample_B.manta_pass.gene_overlap.tsv
    ├── Sample_A.cnv_nondiploid.gene_overlap.tsv
    └── Sample_B.cnv_nondiploid.gene_overlap.tsv
```

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics — 2026-07-20*

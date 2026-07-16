# Project 16 — Wenliang Pan · Human Germline WGS · Analysis Plan

- **Plan Date（创建，immutable）**: 2026-07-15
- **Client**: Wenliang Pan
- **Analyst**: Zhen Gao, PhD — Athenomics
- **Mode**: `/wgs` Mode A — standard germline WGS (nf-core/sarek)
- **Source data**: `/home/gao/Dropbox/Quote_06202601_Wenliang_Pan/`

### 更新记录 / change-log
- 2026-07-15 — 初稿：项目创建、样本表、模式判定、sarek 参数、注释/稀有变异/HLA/报告步骤。

---

## 1. Sample information & input data volume

Canonical sample table: [`sample_info.tsv`](sample_info.tsv) (single source of truth).
Sizes measured from disk `ls -l` on **2026-07-15**; source path above.

| sample | client | species | seq | machine | R1 (GiB) | R2 (GiB) | per-sample (GiB) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Sample_A | Wenliang Pan | Human | WGS PE150 | NovaSeq X Plus | 17.58 | 17.67 | 35.26 |
| Sample_B | Wenliang Pan | Human | WGS PE150 | NovaSeq X Plus | 21.56 | 21.59 | 43.15 |

- **Total dataset: 2 samples, gzip FASTQ 合计 78.40 GiB (`du -shc` = 79 G).**
- Single flowcell/lane `CKDL260012882-1A_23JCJ2LT3_L4`. Estimated depth ~30–40× (confirm from mosdepth after alignment).
- **Provenance check (done):** both files' MD5s match `01.RawData/Sample_A|Sample_B/` in the delivery
  manifest `MD5.txt`; the many `HFD_*/ISRIB_*/Tg_*` entries in that manifest are **other customers'
  samples multiplexed on the same lane**, not ours. No mislabeling. Sample_A/B are our only two.
- **Species sanity check queued:** the lane is dominated by what look like mouse metabolic samples,
  so although the client states Human, we will confirm **GRCh38 primary mapping rate (expect >95%)**
  from sarek/samtools flagstat before trusting human-specific (ClinVar/gnomAD/HLA) annotation.

## 2. Objectives (from client email)

Standard germline WGS: (1) raw-data QC, (2) alignment to human reference, (3) SNV+indel calling,
(4) SV + CNV analysis, (5) annotation vs ClinVar/gnomAD/other DBs, (6) rare + potentially functional
variant identification, (7) HLA typing if feasible, (8) summary report of variants of biological interest.

## 3. Pipeline & step map

| # | step | tool | script |
| :--- | :--- | :---: | :--- |
| 0 | download gnomAD-AF + ClinVar DBs | wget + bcftools | `scripts/0_prep_annotation_dbs.sh` |
| 1 | build sarek samplesheet | python | `scripts/1_make_samplesheet.py` → `scripts/samplesheet.csv` |
| 2 | QC + align + dedup + SNV/indel + SV + CNV + VEP | **nf-core/sarek 3.8.1** | `scripts/2_run_sarek.sh` |
| 3 | run status check | tail | `scripts/3_monitor.sh` |
| 4 | add gnomAD AF + ClinVar to calls | bcftools annotate | `scripts/4_annotate_gnomad_clinvar.sh` |
| 5 | rare + functional variant prioritisation | python | `scripts/5_rare_functional_filter.py` |
| 6 | HLA typing (if feasible) | T1K | `scripts/6_hla_typing.sh` |
| 7 | client report | R/py + markdown | (written at delivery) |

**sarek covers items 1–5 of the objectives in one run:**
- **QC**: FastQC + fastp adapter/quality trimming (`--trim_fastq`) + MultiQC aggregate.
- **Alignment**: bwa-mem2 → GATK MarkDuplicates → CRAM, to **GATK.GRCh38** (analysis-set, the reference
  ClinVar/gnomAD/VEP are all keyed to).
- **SNV/indel**: GATK **HaplotypeCaller** (germline standard, per-sample GVCF → genotyped VCF).
- **SV**: **Manta** + **TIDDIT** (two orthogonal callers; TIDDIT is the safety net if Manta stalls).
- **CNV**: **CNVkit** (germline, flat reference — no matched normal).
- **Functional annotation**: **VEP** (consequence, gene, SIFT, PolyPhen, existing dbSNP IDs); cache
  auto-downloaded via `--download_cache`.

## 4. Analysis rationale & decision criteria

- **`--genome GATK.GRCh38`** (not local GENCODE fasta): the analysis-set naming + provided known-sites
  and intervals are what every clinical DB (ClinVar, gnomAD, VEP) is built against; avoids contig-name
  mismatches. iGenomes bundle auto-fetched once.
- **`--skip_tools baserecalibrator`**: NovaSeq X Plus base qualities are already well-calibrated; BQSR
  gain is negligible and it costs ~20% runtime (proj4-confirmed). No project-specific known-sites needed.
- **`--aligner bwa-mem2`**: faster than bwa; RAM covered (60 GB/task).
- **Annotation split (sarek VEP now, gnomAD+ClinVar as post-step 4)**: VEP inside sarek gives
  consequence/impact; gnomAD **AF-only** (GATK resource, ~3 GB) + ClinVar (~0.2 GB) added by
  `bcftools annotate`. Using AF-only gnomAD (vs the >1 TB full sites VCFs) keeps disk trivial while
  giving exactly the allele-frequency field rarity filtering needs.
- **Rare + functional definition (step 5):**
  - **RARE** = gnomAD_AF < 0.001 **or** absent from gnomAD.
  - **FUNCTIONAL** = VEP HIGH/MODERATE consequence (LoF: stop-gain/frameshift/splice-donor·acceptor/
    start-loss; missense; inframe indels) **or** ClinVar Pathogenic/Likely_pathogenic.
  - **Flagged** = RARE **and** FUNCTIONAL. ClinVar P/LP surfaced regardless of frequency (actionable).
  - Only FILTER=PASS variants evaluated.
- **HLA "if feasible"**: T1K on MHC-region (chr6:28–34 Mb) + unmapped reads extracted from the CRAM —
  class I + II from standard WGS depth. Marked feasible pending depth confirmation.

## 5. Resource allocation

- `scripts/local_resources.config`: `queueSize=2`, `cpus=16`/task, `BWAMEM2 memory=60.GB`,
  defensive `CNNSCOREVARIANTS cpus=9`. CLI `--max_memory 120.GB --max_cpus 56`.
- 2 samples × 16 = **32 threads sustained** (≤ 56 cap ✓). Worst-case RAM 2×60 = 120 GB of 125 GB
  (+65 GB swap) — the proj4-proven envelope.
- Disk: `/home` 696 G free hosts `work/` + `output_results/` + VEP cache (WGS work dir ~200–400 GB
  for 2 samples — fine). Shared annotation DBs on `/Work_bio/.../annotation/`.
- Runs **solo** (no co-running pipeline as of 2026-07-15).

## 6. Run order

```bash
conda run -n regular_bioinfo python scripts/1_make_samplesheet.py   # samplesheet.csv
bash scripts/0_prep_annotation_dbs.sh    # (tmux, parallel) download DBs while sarek runs
bash scripts/2_run_sarek.sh              # (auto tmux 'pan_wgs', auto-resume)
# after sarek completes:
bash scripts/4_annotate_gnomad_clinvar.sh
conda run -n regular_bioinfo python scripts/5_rare_functional_filter.py
bash scripts/6_hla_typing.sh
```

## 7. Estimated timeline

~10–20 h for 2 samples at queueSize=2 (bwa-mem2 align + MarkDuplicates ~3–4 h/sample dominate;
HaplotypeCaller scatter-gather; Manta is the wildcard — TIDDIT backs it up). Then <1 h annotation +
filtering; HLA ~0.5 h.

## 8. Known limitations / open items for client confirmation

- **SNV/indel caller** = HaplotypeCaller only (GATK standard). DeepVariant available as an add/second
  caller if the client wants orthogonal confirmation (CPU-only → much slower).
- **CNV without matched normal**: germline CNVkit uses a flat reference; large/mosaic CNVs are reliable,
  small focal events less so. No panel-of-normals.
- **gnomAD** = AF-only resource (v2-derived, genome+exome joint AF). Adequate for rarity; not
  subpopulation-resolved. Full gnomAD v4 per-population AF can be added on request (disk permitting).
- **HLA** flagged feasible; final call after depth confirmed. T1K gives 2–field typing; if higher
  resolution/validation is needed, HLA-LA (graph-based, heavier) is the alternative.
- **Report emphasis**: to confirm with client — general variant catalogue vs a specific
  gene/phenotype focus (the email says "variants of potential biological significance" generically).

## 9. Pre-delivery self-audit checklist (mandatory)

- MultiQC: per-sample mapping rate (>95% human), duplication, coverage (mosdepth) — any outlier ⇒ investigate.
- Confirm GRCh38 mapping rate before trusting human annotation (species sanity check, §1).
- Ti/Tv ratio of SNVs ~2.0–2.1 genome-wide (sanity of SNV calls).
- Spot-check a ClinVar P/LP hit in IGV; confirm gnomAD_AF annotation actually populated.
- Number sanity across the two samples — wildly different variant counts ⇒ investigate before writing up.

---
*Reference implementations: `4_wgs_human_immu/` (Mode A sarek), `/wgs` skill. Signature on delivery:
Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics.*

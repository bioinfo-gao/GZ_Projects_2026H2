# Whole-Genome Sequencing Analysis — Germline Variant, Structural, HLA and Sample-Origin Profiling

**Report Date:** 2026-07-17
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Homo sapiens* (GRCh38 / GATK.GRCh38 build)
**Tissue / Cell type:** not specified by client
**Project:** Wenliang Pan — human germline WGS (2 samples)
**Samples:** Sample_A, Sample_B (paired-end 150 bp, NovaSeq X Plus)

---

## 1. Objectives

Characterise two human whole-genome sequencing (WGS) samples end-to-end:

1. Confirm species and assess raw-data and alignment quality.
2. Call germline short variants (SNV/indel) genome-wide and annotate them for population rarity (gnomAD) and clinical significance (ClinVar).
3. Prioritise a shortlist of rare and potentially functional / clinically relevant variants.
4. Call structural variants (SV) and copy-number variants (CNV).
5. Perform HLA genotyping (class I and class II) from the WGS data.
6. Infer, from the data alone, whether the material is more consistent with primary/germline tissue or a passaged cell line (the client did not specify tissue of origin).

## 2. Key Findings

- **Both samples are human, high quality.** ≥99.95 % of reads map to GRCh38; post-trim Q30 is 97.3–97.4 %. Mean genome coverage is **23.1× (Sample_A)** and **28.1× (Sample_B)** — adequate for confident germline calling.
- **Genome-wide germline variant catalogues were produced**: **5.88 M** PASS SNV/indel for Sample_A and **6.03 M** for Sample_B, each annotated with gene/consequence (VEP), population frequency (gnomAD) and clinical significance (ClinVar).
- **A prioritised shortlist of rare + functional / clinically flagged variants** was generated: **4,294** for Sample_A and **4,342** for Sample_B (of which **235 / 248** carry a ClinVar Pathogenic/Likely-pathogenic annotation). These are screening candidates, not diagnostic calls.
- **Structural and copy-number callsets** were delivered (Manta, TIDDIT, CNVkit). Both genomes are **predominantly diploid**.
- **HLA genotypes** were typed for both samples, but call confidence varies by locus: class II (**HLA-DRB1/DQB1/DPA1/DPB1/DRB3**) and **HLA-C** are well-supported, whereas the classical class I **HLA-A/B** calls are low-confidence at this WGS depth and should be confirmed by a targeted method before use.
- **Sample-origin inference: both samples are most consistent with primary / germline material, NOT a high-aneuploidy cell line.** The genome is ~90–95 % diploid and long runs of homozygosity are low (4.5–5.0 % of the autosome).

## 3. Sample Information

Canonical sample sheet: `sample_info.tsv`. Raw input sizes measured from disk on 2026-07-15 (source: `/home/gao/Dropbox/Quote_06202601_Wenliang_Pan/`).

| sample   |    client    | species | seq type |    machine    | flowcell/lane | R1 (GiB) | R2 (GiB) | per-sample (GiB) |
| :------- | :----------: | :-----: | :-------: | :------------: | :-----------: | :------: | :------: | :--------------: |
| Sample_A | Wenliang Pan |  Human  | WGS PE150 | NovaSeq X Plus | 23JCJ2LT3_L4 |  17.58  |  17.67  |      35.25      |
| Sample_B | Wenliang Pan |  Human  | WGS PE150 | NovaSeq X Plus | 23JCJ2LT3_L4 |  21.56  |  21.59  |      43.15      |

**Total dataset:** 2 samples, gzip FASTQ **≈ 78.4 GiB** combined.

## 4. Analysis Rationale and Decision Criteria

| Step               | Rationale                                                                                                             | Decision criteria / thresholds                                                           |
| :----------------- | :-------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------- |
| Reference          | GATK.GRCh38 (`Homo_sapiens_assembly38.fasta`) is the standard for GATK germline best-practices                      | chosen for HaplotypeCaller compatibility and<br /> downstream annotation resources       |
| Alignment          | bwa-mem2 (fast, exact drop-in for bwa-mem)                                                                            | mapping rate + properly-paired % as QC gates                                             |
| Duplicate marking  | GATK4 MarkDuplicates (metrics emitted inline)                                                                         | duplicates flagged, not removed                                                          |
| Germline SNV/indel | GATK HaplotypeCaller (best-practices)                                                                                 | report FILTER=PASS variants                                                              |
| Rarity             | a variant is**rare** if gnomAD_AF < 0.001 **or** absent from gnomAD                                       | population-frequency screen for candidate variants                                       |
| Functional         | **HIGH/MODERATE** VEP consequence (LoF, missense, splice, …) **or** ClinVar Pathogenic/Likely-pathogenic | keep = ClinVar P/LP,**or** (rare **and** functional)                         |
| SV / CNV           | Manta + TIDDIT (SV), CNVkit (CNV)                                                                                     | complementary callers; germline (no matched normal)                                      |
| HLA                | T1K`hla-wgs` preset from MHC-region + unmapped reads                                                                | class I + II genotypes                                                                   |
| Origin             | data-driven inference (CNV aneuploidy burden + ROH/LOH) because tissue was not specified                              | aneuploidy fraction > 0.15 and/or ROH fraction > 0.10 would suggest a passaged cell line |

## 5. Methods

- **Pipeline:** nf-core/sarek 3.8.1 (`-profile singularity`), aligner `bwa-mem2`, `--genome GATK.GRCh38`, `--trim_fastq`, tools `haplotypecaller,manta,tiddit,cnvkit,vep`.
- **QC / trimming:** fastp; FastQC; MultiQC.
- **Alignment / dedup:** bwa-mem2 → GATK4 MarkDuplicates; coverage by mosdepth; stats by samtools.
- **Germline calling:** GATK4 HaplotypeCaller (scatter-gather over calling intervals).
- **Annotation:** Ensembl VEP (consequence/gene/SIFT/PolyPhen); then bcftools annotate to add gnomAD allele frequency (`af-only-gnomad.hg38`) and ClinVar `CLNSIG`/`CLNDN` (GRCh38 release).
- **Prioritisation:** in-house filter (`5_rare_functional_filter.py`) applying the rarity + functional criteria in §4.
- **SV/CNV:** Manta (diploid SV), TIDDIT, CNVkit (germline mode).
- **HLA:** T1K 1.0.6, `--preset hla-wgs`, IPD-IMGT/HLA reference.
- **Origin inference:** CNVkit `.call.cns` aneuploidy fraction (genome length in segments with integer copy number ≠ 2) and bcftools `roh` autosomal ROH fraction.

## 6. Results

### 6.1 Raw-data and alignment quality

| metric                         |       Sample_A       |       Sample_B       |
| :----------------------------- | :-------------------: | :-------------------: |
| reads before trim              |      594,077,072      |      731,479,202      |
| reads after trim               |      582,683,738      |      716,718,540      |
| Q30 after trim                 |        97.34 %        |        97.44 %        |
| reads mapped                   | 582,408,655 (99.95 %) | 716,409,454 (99.96 %) |
| properly paired                |        97.7 %        |        97.8 %        |
| duplicate rate                 |        14.4 %        |        15.4 %        |
| mismatch error rate            |        0.437 %        |        0.419 %        |
| mean insert size               |       270.9 bp       |       265.6 bp       |
| **mean genome coverage** |   **23.13×**   |   **28.05×**   |

The ≥99.95 % mapping rate to GRCh38 **confirms human origin** and shows no evidence of foreign/contaminant sequence at the whole-genome level.

### 6.2 Germline SNV/indel (FILTER = PASS)

| metric             | Sample_A | Sample_B |
| :----------------- | :-------: | :-------: |
| total PASS records | 5,876,646 | 6,030,881 |
| SNPs               | 4,781,178 | 4,889,447 |
| indels             | 1,101,392 | 1,148,205 |
| multiallelic sites |  124,950  |  141,910  |

These totals are within the expected range for a human genome at this depth.

### 6.3 Annotation and prioritised variants

Every PASS variant carries VEP consequence + gnomAD_AF + ClinVar CLNSIG/CLNDN in the delivered VCFs. Applying the rarity + functional criteria:

| metric                                   |    Sample_A    |    Sample_B    |
| :--------------------------------------- | :-------------: | :-------------: |
| PASS variants evaluated                  |    5,900,861    |    6,058,771    |
| **prioritised (shortlist)**        | **4,294** | **4,342** |
| — ClinVar Pathogenic/Likely-pathogenic  |       235       |       248       |
| — rare + functional (VEP HIGH/MODERATE) |      4,059      |      4,094      |

The two samples give closely comparable shortlist sizes (4,294 vs 4,342), as expected for two human genomes processed identically. The prioritised tables (`*.prioritised.tsv`) list, per variant: `chrom, pos, ref, alt, gnomAD_AF, gene, consequence, clinvar_sig, clinvar_disease, flag`. **These are screening candidates for expert review, not clinical diagnoses**; many ClinVar entries carry "Conflicting classifications" and require manual curation and, where relevant, orthogonal confirmation.

### 6.4 Structural and copy-number variants

| caller                   | Sample_A | Sample_B |
| :----------------------- | :------: | :------: |
| Manta diploid SV (PASS)  |  6,349  |  7,209  |
| TIDDIT SV (all)          |  55,326  |  62,469  |
| CNVkit called segments   |   343   |   303   |
| — segments with CN ≠ 2 |   192   |   162   |

CNVkit was run in germline mode (no matched normal), so copy-number is relative to a flat/pooled baseline; small aberrant segments in low-mappability/segmental-duplication regions are largely technical. Although about half of the called *segments* carry CN ≠ 2, these are mostly short, so by genome *length* only ~5–10 % is non-diploid (see §6.6) — i.e. the genome is predominantly diploid. Per-sample CNV scatter plots are provided (`structural_cnv/*-scatter.png`).

### 6.5 HLA genotyping (T1K, `hla-wgs`)

Both samples were typed at class I and class II loci (37 of 41 T1K loci returned at least one allele). **Call confidence differs markedly by locus and must be read alongside the T1K quality score** (higher = better supported; in this run reliable calls score ≥13, weak calls score 0–1). The table below gives the representative genotype per locus with each allele's T1K quality in parentheses:

| locus    | Sample_A (allele / T1K quality)          | Sample_B (allele / T1K quality)          | confidence  |
| :------- | :--------------------------------------- | :--------------------------------------- | :---------: |
| HLA-A    | A\*23:144 (1) / A\*01:01:01 (0)          | A\*30:02:01 (22) — single                | A: mixed; B: **low** |
| HLA-B    | B\*15:16:01 (0) / B\*35:03:01 (0)        | B\*35:03:01 (1) / B\*58:149 (1)          | **low**     |
| HLA-C    | C\*07:18:01 (21) — single                | C\*07:18:01 (27) — single                | high        |
| HLA-DRB1 | DRB1\*13:03:01 (13) / DRB1\*13:01:01 (4) | DRB1\*13:04 (10) / DRB1\*13:01:01 (5)    | high / mod  |
| HLA-DRB3 | DRB3\*02:02:01 (31) — single             | DRB3\*02:02:01 (34) — single             | high        |
| HLA-DQB1 | DQB1\*03:03:02 (18) / DQB1\*02:02 (15)   | DQB1\*03:19:01 (18) / DQB1\*03:03:02 (18)| high        |
| HLA-DPA1 | DPA1\*02:01:01 (36) / DPA1\*01:03:01 (34)| DPA1\*01:03:01 (33) / DPA1\*02:01:01 (25)| high        |
| HLA-DPB1 | DPB1\*17:01:01 (27) / DPB1\*04:01:01 (20)| DPB1\*11:01:01 (32) / DPB1\*04:01:01 (33)| high        |

Reading the results:

- **Well-supported (use with confidence):** the class II loci **DRB1, DQB1, DPA1, DPB1, DRB3** (quality 13–36) and **HLA-C** (21–27), plus **HLA-A in Sample_B** (22). The class II calls are internally coherent — both samples carry **DRB1\*13 + DRB3\*02:02 and no DRB5**, the canonical DR13 / DR52 haplotype linkage — which corroborates their validity.
- **Low-confidence (confirm before use):** the classical class I **HLA-A (Sample_A) and HLA-B (both samples)** scored quality **0–1** with low read support (abundance ~1–4 vs ~30–50 for the reliable calls). This is the expected limitation of WGS at this depth over the hyper-polymorphic class I exons 2/3; the **specific allele digits are not reliable** and should be confirmed with a targeted method (amplicon- or long-read-based HLA typing) if clinically required.
- **Single-allele loci** (HLA-C in both samples; HLA-A in Sample_B) may be true homozygotes **or** reflect second-allele dropout — the two cannot be distinguished from this data.
- **Allele ambiguity:** **DRB4** and some second alleles of **DQA1 / DQB1** are reported as comma-separated candidate lists in the genotype file (the gene is present but the exact allele is unresolvable); the first listed allele is the representative.

Full per-locus genotypes with quality and abundance are in `hla_typing/*_hla_genotype.tsv`.

### 6.6 Sample-origin inference

| signal                                       |   Sample_A   |   Sample_B   | interpretation                          |
| :------------------------------------------- | :----------: | :----------: | :-------------------------------------- |
| aneuploidy fraction (genome length, CN ≠ 2) |    0.104    |    0.049    | predominantly diploid                   |
| autosomal ROH fraction                       |    0.050    |    0.045    | low LOH                                 |
| V(D)J clonotypes                             | not assessed | not assessed | TRUST4 reference unavailable this round |

Both samples are **most consistent with primary / germline material**. The genome is ~90–95 % diploid and ROH is low — this is **not** the extensive aneuploidy or long LOH tracts characteristic of a passaged/clonal cell line. This is a data-driven inference and **does not replace the client stating the actual tissue/cell of origin**; V(D)J-based lymphoid profiling was not assessable this round (reference not configured) and can be added on request.

## 7. Conclusions

| # | Conclusion                                                             | Evidence                                                                  |
| :- | :--------------------------------------------------------------------- | :------------------------------------------------------------------------ |
| 1 | Both samples are human and of high quality                             | ≥99.95 % mapping to GRCh38; Q30 97.3–97.4 %; 23–28× coverage          |
| 2 | Complete, annotated germline SNV/indel catalogues delivered            | 5.88 M / 6.03 M PASS variants with VEP + gnomAD + ClinVar                 |
| 3 | Prioritised rare/functional shortlists delivered for expert review     | 4,294 / 4,342 variants (235 / 248 ClinVar P/LP)                           |
| 4 | SV and CNV callsets delivered; genomes predominantly diploid           | Manta/TIDDIT/CNVkit; CN ≠ 2 over ~5–10 % of genome length               |
| 5 | HLA typed; class II + HLA-C well-supported, class I HLA-A/B low-confidence | T1K`hla-wgs`; per-locus quality tiering in §6.5                     |
| 6 | Origin most consistent with primary/germline material, not a cell line | aneuploidy 0.05–0.10, ROH 0.045–0.050 (both below cell-line thresholds) |

**Caveats.** Prioritised variants and ClinVar annotations are screening-level and require expert curation; germline CNV (no matched normal) is baseline-relative; HLA class I (especially HLA-A/B) is low-confidence at this depth and needs orthogonal confirmation, and some loci carry allele ambiguity; origin is an inference, not a substitute for known provenance.

## 8. Deliverable Files

```
custom_research_report_20260717/
├── Wenliang_Pan_WGS_0717.md            ← this report
├── qc/                                  MultiQC + fastp HTML
├── alignment_coverage/                  mosdepth summaries + samtools stats
├── variant_calling_snv_indel/           HaplotypeCaller PASS VCFs + VEP-annotated VCFs (+ .tbi)
├── annotation_prioritised/              gnomAD+ClinVar VCFs; per-sample & combined prioritised .tsv
├── structural_cnv/                      Manta + TIDDIT VCFs, CNVkit .call.cns + scatter PNGs
├── hla_typing/                          per-sample HLA genotype .tsv
└── sample_origin/                       origin_summary.tsv + per-sample ROH
```

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics.*

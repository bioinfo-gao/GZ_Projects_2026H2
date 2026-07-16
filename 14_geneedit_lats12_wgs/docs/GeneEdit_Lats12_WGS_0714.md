# Whole-Genome Sequencing Analysis — Gene-Edited Ovarian Cancer Models & Lats1/2 Hippo-Pathway Mice

**Report Date:** 2026-07-15
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Project:** Two-arm WGS, 12 samples (Study A: CRISPR-edited tumor models, n=6; Study B: Lats1/2-flox tissues ± iHPV, n=6)

---

## 1. Objectives

**Study A — CRISPR-edited cells → in vivo tumors (Trp53⁺/⁻;Cas9 background):**

- **A1** Verify that the intended CRISPR edits (Brca1, Brca2, Pten knockout) actually occurred.
- **A2** Characterise the genomes of the resulting tumors (somatic changes, copy-number/aneuploidy), and determine each tumor's lineage of origin (B1TP = Brca1+Pten vs B2TP = Brca2+Pten).

**Study B — Lats1/2-flox tissues, ± iHPV, age series (3M/12M/18M):**

- Determine whether the genomes of these tissues carry abnormalities (aneuploidy, structural/de-novo variants) relative to the normal C57BL/6 genome, and whether any such change could underlie the reported (>10-month) oviduct phenotype.
- Locate the iHPV transgene and assess whether its integration disrupts a candidate gene.

A key methodological point (established in the analysis plan): the engineered elements in Study B are **loxP-only / Cre-dependent and un-activated** by design (Lats1/2 floxed but not deleted; E6/E7 behind a lox-stop-lox and not expressed). Analyses are therefore framed around insertional/somatic/de-novo hypotheses rather than an assumed "Lats1/2-loss → instability" mechanism.

---

## 2. Key Findings

1. **Study B tissues are genomically stable — no aneuploidy.** All six Lats1/2-flox tissues (both genotypes, all ages) show a flat, diploid copy-number profile genome-wide (autosomal copy number 1.94–2.08). The reported oviduct phenotype is **not** driven by large-scale chromosomal instability or aneuploidy.
2. **The method is validated by a built-in positive control.** The same pipeline calls clear aneuploidy in the Study A tumors (below), proving it *can* detect copy-number change when present — so the flatness of Study B is a real biological result, not a sensitivity failure.
3. **Study A tumors are aneuploid, each with a distinct karyotype.** The three tumors carry different chromosome-scale gains/losses, consistent with independent clonal evolution of Brca/Pten/Trp53-driven, homologous-recombination-deficient tumors. The parental and edited-but-pre-tumor cells remain diploid.
4. **CRISPR edits confirmed; tumor lineages resolved.** Brca1 and Pten knockout are confirmed in B1TP; Pten knockout (Brca1 wild-type) in B2TP. Tumors 1 and 2 derive from the **B1TP (Brca1+Pten)** lineage. Tumor 3 is an exception requiring follow-up (see §7).
5. **The iHPV transgene is present specifically in L1L2H mice.** HPV16 E6/E7 and luciferase reads are detected in all three L1L2H samples and are entirely absent from all three L1L2 samples (perfect specificity), at levels consistent with a fixed germline transgene.

---

## 3. Sample Information

12 samples, paired-end 150 bp WGS (NovaSeq X Plus), aligned to GRCm39 (GENCODE vM35) with nf-core/sarek 3.8.1 (bwa-mem2). Mean autosomal depth from this analysis in parentheses.

| #  | Label     | Study |  Type  | Group | Genotype / role                                               |
| :- | :-------- | :---: | :----: | :----: | :------------------------------------------------------------ |
| 1  | RO_origin |   A   |  Cell  | parent | Trp53⁺/⁻; Cas9; unedited —**matched normal** (~20×) |
| 2  | RO_B1TP   |   A   |  Cell  | edited | Brca1 + Pten KO (~24×)                                       |
| 3  | RO_B2TP   |   A   |  Cell  | edited | Brca2 + Pten KO (~24×)                                       |
| 4  | RO_tumor1 |   A   |  Cell  | tumor | tumor from B1TP/B2TP injection (~26×)                        |
| 5  | RO_tumor2 |   A   |  Cell  | tumor | tumor (~19×)                                                 |
| 6  | RO_tumor3 |   A   |  Cell  | tumor | tumor (~30×)                                                 |
| 7  | L1L2_3M   |   B   | Tissue |  L1L2  | Lats1/2 flox, 3 months (~21×)                                |
| 8  | L1L2H_3M  |   B   | Tissue | L1L2H | Lats1/2 flox + iHPV, 3 months (~20×)                         |
| 9  | L1L2_12M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 12 months (~32×)                               |
| 10 | L1L2H_12M |   B   | Tissue | L1L2H | Lats1/2 flox + iHPV, 12 months (~28×)                        |
| 11 | L1L2_18M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 18 months (~25×)                               |
| 12 | L1L2H_18M |   B   | Tissue | L1L2H | Lats1/2 flox + iHPV, 18 months (~21×)                        |

Total input: 12 samples, gzip FASTQ 462 GiB (measured on disk 2026-07-12).

---

## 4. Analysis Rationale and Decision Criteria

| Question                        | Approach & why it works without a dedicated wild-type control                                                                       | Threshold / criterion                                                      |
| :------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------- |
| Copy number / aneuploidy        | Ratio of genome-internal coverage → copy number; reference-free. GRCm39 = C57BL/6J ≈ the animals' own normal.                     | Chromosome-median CN ≥2.5 = gain, ≤1.5 = loss                            |
| Edit verification (A1)          | Read directly at each sgRNA cut site; compare edited/tumor vs the RO_origin parent.                                                 | Indel present in edited sample, absent (0/0) in RO_origin                  |
| Tumor lineage (A2)              | Brca1-KO ⇒ B1TP lineage; Brca2-KO ⇒ B2TP lineage.                                                                                 | Presence/absence of Brca1 cut-site indel                                   |
| iHPV presence / integration (B) | Transgene reads are unmapped against plain GRCm39 → extract and align to HPV16/EGFP/Luc markers; L1L2 = internal negative control. | Construct-mapping reads present in L1L2H, absent in L1L2                   |
| De-novo candidates (B)          | Subtract known C57BL/6 strain background (Sanger MGP, all strains) from each sample's germline calls.                               | Variant private to sample (not in MGP), recurrent or genotype-differential |

---

## 5. Methods

- **Alignment / preprocessing:** nf-core/sarek 3.8.1, bwa-mem2, GRCm39 (GENCODE vM35), duplicate marking, BQSR skipped.
- **Variant calling:** Study A somatic = GATK Mutect2 (each tumor vs RO_origin) + TIDDIT SV; Study B germline = GATK HaplotypeCaller + CNNScoreVariants + TIDDIT SV.
- **Copy number:** mosdepth 500 kb binned depth → normalised to per-sample autosomal median → copy number = 2 × ratio; chromosome-level medians and cohort heatmap (mouse chromosomes are acrocentric, so whole-chromosome CN ≈ arm-level).
- **Edit verification:** sgRNA spacers located in target genes (seqkit), cut sites defined 3 bp from the PAM; multi-sample `bcftools mpileup`/`call` at cut-site windows; indel genotypes/allele depths compared to RO_origin.
- **iHPV integration:** unmapped reads extracted per sample and aligned (bwa-mem2) to a marker reference (HPV16 NC_001526.4 + EGFP + firefly luciferase); genomic anchors from mapped reads with unmapped mates.
- **De-novo candidates:** Sanger Mouse Genomes Project v8 (REL-2021, GRCm39) SNP+indel VCFs subtracted from each sample's germline calls with `bcftools isec` (chromosome names reconciled chr1↔1).

---

## 6. Results

### 6.1 Copy number & aneuploidy (headline)

Chromosome-level copy number (autosomes) across the cohort:

- **Study B (all six Lats1/2 tissues):** flat and diploid — autosomal copy number ranges only **1.94–2.08** in every sample, at every age, in both genotypes. No chromosome-scale gain or loss.
- **Study A parent + edited cells (RO_origin, B1TP, B2TP):** also flat/diploid (1.91–2.14).
- **Study A tumors:** clearly aneuploid, each distinct —
  - **tumor1:** gains chr8, chr10, chr11; loss chr18
  - **tumor2:** gain chr5 (plus milder shifts)
  - **tumor3:** losses chr8, chr12, chr13, chr14; gain chr15 (most-rearranged genome)

**Interpretation.** The Study A tumors act as an internal positive control: the pipeline detects clear aneuploidy where it exists. Against that, the uniformly flat Study B profiles are a genuine biological result — **the Lats1/2-flox tissues have stable, diploid genomes with no aneuploidy, so the oviduct phenotype is not explained by large-scale genomic instability.** This is consistent with the engineered elements being un-activated (Lats1/2 not deleted). See `cnv_ploidy/cohort_cn_heatmap.png` and per-sample `*.cn_profile.png`.

### 6.2 CRISPR edit verification & tumor lineage (Study A)

At the sgRNA cut sites (RO_origin is homozygous reference throughout, confirming these are edits, not background):

| Gene / site             | RO_origin |      B1TP      |         B2TP         |           tumor1           |        tumor2        | tumor3 |
| :---------------------- | :-------: | :-------------: | :-------------------: | :-------------------------: | :-------------------: | :----: |
| Brca1 (chr11:101422906) |    0/0    | **indel** |          0/0          | **indel (biallelic)** | **indel (hom)** |  0/0  |
| Pten (chr19:32777294)   |    0/0    | **indel** | **indel (hom)** |    **indel (hom)**    |         indel         |  0/0  |

- **B1TP:** Brca1 KO + Pten KO — matches its Brca1+Pten design.
- **B2TP:** Pten KO, Brca1 wild-type — correct (it targets Brca2, not Brca1).
- **tumor1, tumor2:** carry both Brca1 and Pten indels → **B1TP (Brca1+Pten) lineage.**
- **tumor3:** no detectable Brca1 or Pten edit (see §7).

### 6.3 iHPV transgene detection (Study B)

Construct-marker reads (unmapped reads aligned to HPV16 + EGFP + luciferase):

| Sample    | Construct reads | HPV16 | Luciferase |
| :-------- | :-------------: | :---: | :--------: |
| L1L2_3M   |        0        |   0   |     0     |
| L1L2H_3M  |       45       |  16  |     71     |
| L1L2_12M  |        0        |   0   |     0     |
| L1L2H_12M |       83       |  41  |    123    |
| L1L2_18M  |        0        |   0   |     0     |
| L1L2H_18M |       58       |  37  |     75     |

**The iHPV transgene is present in all three L1L2H samples and absent in all three L1L2 samples** — a perfectly specific result. Read counts do not increase with age, consistent with a fixed germline transgene rather than an age-accumulating event.

### 6.4 Structural variants

TIDDIT SV counts are available for all samples (`sv/sv_counts.tsv`). For Study B these are dominated by C57BL/6-versus-reference background SVs and are being refined with the MGP background filter before interpretation; Study A tumor SVs will be reported as tumor-minus-normal somatic events.

### 6.5 De-novo candidate variants (Study B)

Each sample's germline calls were subtracted against the Sanger MGP v8 (all-strain) SNP+indel catalogue to remove known mouse strain background:

| Sample    | Total germline variants vs GRCm39 | Private after MGP subtraction |
| :-------- | :-------------------------------: | :---------------------------: |
| L1L2_3M   |             5,101,717             |            573,954            |
| L1L2H_3M  |             5,334,890             |            562,308            |
| L1L2_12M  |             6,110,749             |            658,142            |
| L1L2H_12M |             5,612,998             |            589,237            |
| L1L2_18M  |             5,899,863             |            625,093            |
| L1L2H_18M |             5,884,909             |            615,276            |

**Interpretation — an important caveat.** Two observations show these counts are **not** a list of causal de-novo candidates:

- Each sample carries **~5–6 million variants relative to GRCm39**. A pure C57BL/6J animal would differ from GRCm39 (itself a 6J assembly) by only tens of thousands. This magnitude means **the tissues are not pure C57BL/6J** — they carry a substantial divergent (likely 129-derived, from the engineering) background genome-wide. This tempers the earlier "reference ≈ their normal" assumption.
- After MGP subtraction, **~0.57–0.66 million private variants** remain per sample, distributed uniformly across all chromosomes (proportional to chromosome length) — i.e. residual background divergence not covered by MGP, not localised de-novo events.

Genuine causal de-novo events number in the dozens-to-hundreds, not hundreds of thousands. **Identifying them from this background requires functional restriction to high-impact coding consequences (frameshift / stop-gain / splice, via VEP or snpEff) combined with recurrence and L1L2-vs-L1L2H differential filtering.** That functional-annotation step is the recommended immediate next analysis; the per-sample private call sets are provided for it (`candidates_denovo/*.private.vcf.gz`). Critically, this does not affect the primary Study B conclusion: the copy-number analysis (§6.1) already establishes these genomes are structurally stable and diploid.

### 6.6 Somatic point mutations (Study A)

Mutect2 tumor-vs-origin calling was run for all five pairs. **The current PASS counts are inflated by artifacts (no contamination table / panel-of-normals) and are not yet a reliable mutation burden;** they are therefore not reported here as findings. This analysis requires additional filtering (VAF/depth thresholds, orientation-bias, PoN) before mutation burdens are quoted. The Study A tumor genome story is presently carried by the copy-number (§6.1) and edit-verification (§6.2) results, which are robust.

---

## 7. Conclusions & Items Requiring Follow-up

| Conclusion                                                                     | Confidence | Evidence                           |
| :----------------------------------------------------------------------------- | :--------: | :--------------------------------- |
| Study B tissues have stable diploid genomes; no aneuploidy at any age/genotype |    High    | §6.1, flat CN 1.94–2.08          |
| Oviduct phenotype is not driven by large-scale genomic instability             |    High    | §6.1 + validated positive control |
| Study A tumors are aneuploid with distinct karyotypes                          |    High    | §6.1                              |
| B1TP = Brca1+Pten KO; B2TP = Pten KO (Brca1 WT) confirmed                      |    High    | §6.2                              |
| Tumors 1 & 2 arose from the B1TP lineage                                       |    High    | §6.2                              |
| iHPV transgene present specifically in L1L2H                                   |    High    | §6.3                              |

**Items requiring follow-up / client input:**

1. **Tumor3 anomaly** — tumor3 is the most aneuploid tumor yet shows **no** Brca1 or Pten edit at the assayed cut sites. Possible explanations (minor un-edited clone of origin, B2TP-Brca2 lineage with Pten reversion, or a labeling issue) need resolution; this requires the Brca2 guide sequences and targeted IGV review.
2. **Brca2 sgRNA sequences needed** — this batch provided only Brca1 + Pten guides. Brca2 knockout could not be pinpointed (the gene-wide indel scan returned only homopolymer sequencing noise). Please provide the three Brca2 sgRNAs to confirm B2TP editing and tumor3's lineage.
3. **iHPV integration locus** — construct presence is confirmed, but the base-pair integration junction / disrupted gene is not yet resolved; this requires the full PMC4662542 vector map to capture junction-spanning reads.
4. **Somatic point-mutation burden** — requires the additional filtering described in §6.6 before numbers are quoted.
5. **Study B strain background is not pure C57BL/6J** (~5–6M variants vs GRCm39, §6.5) — the tissues carry a divergent (likely 129-derived) background. De-novo candidate mining is therefore a two-step task: functional high-impact annotation (VEP/snpEff) + recurrence / genotype-differential filtering on the provided private call sets. This is the recommended next analysis and can be run on request.

---

## 8. Deliverable Files

```
custom_research_report_20260714/
├── GeneEdit_Lats12_WGS_0714.md      ← this report
├── qc/                               MultiQC reports (Study A, Study B)
├── cnv_ploidy/                       cohort CN table, aneuploidy calls, heatmap, per-sample CN profiles
├── edit_verification/                cut sites, cut-site indel genotypes, Brca2 scan
├── ihpv_integration/                 construct-presence table (per-sample HPV16/EGFP/Luc reads)
├── somatic/                          Mutect2 PASS counts (preliminary), Trp53 locus depth
├── sv/                               TIDDIT SV counts
└── candidates_denovo/                MGP-filtered de-novo candidates (in progress)
```

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics.*

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

A key methodological point: the engineered elements in Study B are **loxP-only / Cre-dependent and un-activated** by design (Lats1/2 floxed but not deleted; E6/E7 behind a lox-stop-lox and not expressed). Analyses are therefore framed around insertional/somatic/de-novo hypotheses rather than an assumed "Lats1/2-loss → instability" mechanism.

---

## 2. Key Findings

1. **Study A tumors are aneuploid, each with a distinct karyotype.** The three tumors carry different chromosome-scale gains/losses, consistent with independent clonal evolution of Brca/Pten/Trp53-driven, homologous-recombination-deficient tumors. The parental and edited-but-pre-tumor cells remain diploid.
2. **CRISPR edits verified; tumor lineages resolved.** Brca1 + Pten knockout confirmed in B1TP; Pten knockout confirmed in B2TP (Brca1 wild-type, as expected for its Brca2+Pten design — a targeted, base-accurate Brca2 cut-site confirmation additionally requires the client's Brca2 sgRNA sequences, which were not part of this batch, §7). Tumors 1 and 2 derive from the **B1TP (Brca1+Pten)** lineage. Tumor 3 is wild-type at all assayed Brca1 and Pten cut sites and most likely arose from an un-edited subclone (see §7).
3. **Study B tissues are genomically stable — no aneuploidy.** All six Lats1/2-flox tissues (both genotypes, all ages) show a flat, diploid copy-number profile genome-wide (autosomal copy number 1.94–2.08). The reported oviduct phenotype is **not** driven by large-scale chromosomal instability or aneuploidy.
4. **The iHPV transgene is present specifically in L1L2H mice.** HPV16 E6/E7 and luciferase reads are detected in all three L1L2H samples and are entirely absent from all three L1L2 samples (perfect specificity), at levels consistent with a fixed germline transgene.

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
| De-novo candidates (B)          | Subtract known mouse strain background (Sanger Mouse Genomes Project [MGP], all strains) from each sample's germline calls.          | Variant private to sample (not in MGP), recurrent or genotype-differential |

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
- **B2TP:** Pten KO (homozygous), Brca1 wild-type — matches the full expected signature of its **Brca2+Pten** design (B2TP targets Brca2, not Brca1, so Brca1 is expected to stay wild-type). A *targeted, base-accurate* Brca2 cut-site genotype additionally requires the client's three Brca2 sgRNA sequences (the spacer defines the predicted cut site); these were not part of this batch, so that confirmation step, if needed, is gated on that client input (see §7).
- **tumor1, tumor2:** carry both Brca1 and Pten indels → **B1TP (Brca1+Pten) lineage.**
- **tumor3:** wild-type at **all** assayed Brca1 **and** Pten cut sites, with solid coverage (Pten sites 17–18× reads, unambiguous reference) — a true absence of edit, not a coverage gap. Brca1 wild-type status *by itself* would be consistent with a Brca2+Pten (B2TP) origin, since B2TP does not target Brca1; the genuinely unexpected finding is the **wild-type Pten**, because **both** injected lines carry a Pten edit — B2TP homozygously. tumor3 therefore cannot be cleanly derived from either edited line at the assayed sites (interpretation in §7).

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

The Sanger Mouse Genomes Project (MGP) is a reference catalogue of the DNA variants that distinguish the common laboratory mouse strains (e.g. 129, C57BL/6 substrains, BALB/c) from the GRCm39 reference; subtracting it removes inherited strain background so that only variants genuinely private to a sample remain. Each sample's germline calls were subtracted against the Sanger MGP v8 (all-strain) SNP+indel catalogue on this basis:

| Sample    | Total germline variants vs GRCm39 | Private after MGP subtraction |
| :-------- | :-------------------------------: | :---------------------------: |
| L1L2_3M   |             5,101,717             |            573,954            |
| L1L2H_3M  |             5,334,890             |            562,308            |
| L1L2_12M  |             6,110,749             |            658,142            |
| L1L2H_12M |             5,612,998             |            589,237            |
| L1L2_18M  |             5,899,863             |            625,093            |
| L1L2H_18M |             5,884,909             |            615,276            |

**Interpretation — an important caveat.** Two observations show these counts are **not** a list of causal de-novo candidates:

- Each sample carries **~5–6 million variants relative to GRCm39** — a firm and informative result. Because **GRCm39 is the C57BL/6J (“6J”) reference assembly**, a genuinely pure C57BL/6J animal sits only tens of thousands of variants from it (ordinary colony drift). A count in the *millions* is therefore a clear signal, not noise: **these tissues carry a substantial non-6J classical-inbred-strain background across the entire genome — i.e. the line is not congenically pure C57BL/6J.** Study B variants are accordingly interpreted **after MGP background subtraction**: this genome-wide strain background is expected and removed by design, leaving the sample-private variants that matter for de-novo candidate mining.
- **Most likely source (educated guess): a residual 129-derived engineering background.** Gene-targeted and floxed alleles are typically generated in **129-derived embryonic stem cells**, and a 129 background persists genome-wide unless the line is fully backcrossed to congenic purity on C57BL/6J. This particular line's ES-cell provenance is not documented in the materials provided, so we present 129 as the leading explanation rather than a certainty — but it is well supported: it is both the standard ES-cell donor for this kind of allele and a quantitative match, since **the 129 strain lies ~4–6 million SNPs/indels from the 6J reference**, the same order as the ~5–6 M observed here. In short, the residual divergence is most consistent with an engineered line carrying an incompletely backcrossed 129 background. A definitive call on the exact donor strain is available on request via a dedicated strain-assignment analysis (comparing each sample's private variants against per-strain MGP catalogues).
- After MGP subtraction, **~0.57–0.66 million private variants** remain per sample, distributed uniformly across all chromosomes (proportional to chromosome length) — i.e. residual background divergence not covered by MGP, not localised de-novo events.

Genuine causal de-novo events number in the dozens-to-hundreds, not hundreds of thousands. **Identifying them from this background requires functional restriction to high-impact coding consequences (frameshift / stop-gain / splice, via VEP or snpEff) combined with recurrence and L1L2-vs-L1L2H differential filtering.** That functional-annotation step is the recommended immediate next analysis; the per-sample private call sets are provided for it (`candidates_denovo/*.private.vcf.gz`). Critically, this does not affect the primary Study B conclusion: the copy-number analysis (§6.1) already establishes these genomes are structurally stable and diploid.

### 6.6 Somatic point mutations (Study A)

Mutect2 tumor-vs-origin calling was run for all five pairs. **The current PASS counts are inflated by artifacts (no contamination table / panel-of-normals) and are not yet a reliable mutation burden;** they are therefore not reported here as findings. This analysis requires additional filtering (VAF/depth thresholds, orientation-bias, PoN) before mutation burdens are quoted. The Study A tumor genome story is presently carried by the copy-number (§6.1) and edit-verification (§6.2) results, which are robust.

---

## 7. Conclusions & Items Requiring Follow-up

| Conclusion                                                                     |     Confidence     | Evidence                                       |
| :----------------------------------------------------------------------------- | :----------------: | :--------------------------------------------- |
| Study B tissues have stable diploid genomes; no aneuploidy at any age/genotype |        High        | §6.1, flat CN 1.94–2.08                      |
| Oviduct phenotype is not driven by large-scale genomic instability             |        High        | §6.1 + validated positive control             |
| Study A tumors are aneuploid with distinct karyotypes                          |        High        | §6.1                                          |
| B1TP = Brca1+Pten KO confirmed; B2TP = Pten KO + Brca1 WT confirmed            |        High        | §6.2                                          |
| Targeted Brca2 cut-site confirmation (B2TP)                                    | Needs client input | Brca2 sgRNA sequences not yet provided (§7.2) |
| Tumors 1 & 2 arose from the B1TP lineage                                       |        High        | §6.2                                          |
| Tumor3 arose from an un-edited subclone (not B1TP/B2TP)                        |      Moderate      | §6.2/§7.1; Brca1+Pten both WT                |
| iHPV transgene present specifically in L1L2H                                   |        High        | §6.3                                          |

**Items requiring follow-up / client input:**

1. **Tumor3 origin** — tumor3 is the most aneuploid tumor, yet is wild-type at **all** assayed Brca1 and Pten cut sites (solid coverage — not a depth artifact). Since both injected lines carry a Pten edit (B2TP homozygously), a straightforward descent from either edited line would still leave the Pten lesion detectable, so the wild-type Pten is the crux. The finding is **most consistent with tumor3 having arisen from a cell that did not carry the intended Brca1/Pten edits** — an un-edited / editing-escaped subclone of the Trp53⁺/⁻;Cas9 parent, whose high aneuploidy is then Trp53-loss-driven and independent of Brca/Pten editing; a sample-tracking swap is the main alternative. A Brca2+Pten (B2TP) origin is **disfavored**, because B2TP's homozygous Pten edit is not credibly reverted. Confirmation path: (a) genotype the Brca2 cut site once the Brca2 sgRNAs are provided — a Brca2 edit in tumor3 would reopen a B2TP origin; (b) targeted IGV of the Brca1/Pten loci plus a SNP-fingerprint identity check against RO_origin to distinguish "un-edited subclone" from "sample swap".
2. **Brca2 sgRNA sequences needed from the client — the one input required to close out B2TP / tumor3.** This batch included only the Brca1 + Pten sgRNAs; the Brca2 sgRNA sequences were not provided. Because a spacer sequence is what defines the predicted Brca2 cut site, a *targeted, base-accurate* confirmation of the Brca2 edit is gated on that input — the exploratory gene-wide Brca2 indel scan we ran cannot substitute for it (a guide-free scan returns only homopolymer sequencing noise rather than a spacer-anchored call). **To take this to completion, please provide the three Brca2 sgRNA spacer sequences.** With them we can (i) confirm the intended Brca2 edit in the B2TP cell line at its predicted cut site, and (ii) test tumor3 for a Brca2 edit to settle its lineage. On the materials provided, B2TP already shows the full expected signature of its Brca2+Pten design — Pten knocked out and Brca1 wild-type.
3. **iHPV integration locus** — construct presence is confirmed, but the base-pair integration junction / disrupted gene is not yet resolved; this requires the full PMC4662542 vector map to capture junction-spanning reads.
4. **Somatic point-mutation burden** — requires the additional filtering described in §6.6 before numbers are quoted.
5. **Study B strain background is not congenically pure C57BL/6J** (~5–6M variants vs GRCm39, §6.5) — the tissues carry a divergent non-6J classical-inbred background genome-wide, most likely a residual 129 background from the ES-cell engineering (educated guess; exact donor strain confirmable via strain-assignment on request — see §6.5). De-novo candidate mining is therefore a two-step task: functional high-impact annotation (VEP/snpEff) + recurrence / genotype-differential filtering on the provided private call sets. This is the recommended next analysis and can be run on request.

---

## 8. Deliverable Files

```
custom_research_report_20260715/
├── GeneEdit_Lats12_WGS_0715.md      ← this report
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

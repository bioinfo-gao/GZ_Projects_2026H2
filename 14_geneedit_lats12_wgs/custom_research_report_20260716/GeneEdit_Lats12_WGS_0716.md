# Whole-Genome Sequencing Analysis — Gene-Edited Ovarian Cancer Models & Lats1/2 Hippo-Pathway Mice

**Report Date:** 2026-07-16
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Mus musculus* (GRCm39, GENCODE vM35)
**Tissue / Cell:** Study A — cultured cells: primary cells from a Trp53⁺/⁻;Cas9 mouse and cell lines digested from the resulting in-vivo tumors. Study B — mouse tissue (oviduct-associated; specific tissue not specified by client)
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

1. **All 12 samples pass QC** — 99.99% of reads map to GRCm39, duplicate rates 7.7–14.3%, mean depth 19.2–31.9× with 84–92% of the genome at ≥10×. No sample was excluded and no conclusion below is limited by data quality (§6.1).
2. **All intended CRISPR edits are verified, and the tumor lineages are resolved.** Brca1 + Pten knockout confirmed in B1TP. **Brca2 knockout confirmed in B2TP, and it is biallelic** — no wild-type Brca2 allele survives — alongside its Pten knockout. (The knockout itself is solid; *which* lesion sits on each Brca2 allele is not resolvable at this depth and is not claimed — §6.3.) Tumors 1 and 2 derive from the **B1TP (Brca1+Pten)** lineage.
3. **Tumor 3 carries none of the intended edits and arose from an un-edited subclone of the parent.** It is wild-type at all three targeted genes (Brca1, Brca2, Pten), yet a SNP-fingerprint check shows it is genetically the **same animal** as the parental line — so it is an editing-escaped subclone, not a mis-tracked sample. **Practical consequence: editing-escaped cells can still form tumors in this model** (tumor3 is in fact the most aneuploid of the three), so a tumor phenotype should not be attributed to Brca/Pten loss without genotyping that specific tumor.
4. **All three tumors lost Trp53 — the second hit that drove them.** The parent is Trp53⁺/⁻; every tumor has become homozygous at Trp53, while the edited pre-tumor cells have not. In tumor1 the homozygous tract is **focal, centred on Trp53** (92% at 69–70 Mb, against 1–3% just 3 Mb either side) — selection for Trp53 loss, visible directly in the data (§6.5).
5. **Tumor3 has undergone a whole-genome catastrophe: 95% of its genome has lost heterozygosity, at diploid copy number.** This is genome-wide copy-neutral LOH — loss of one haploid genome followed by duplication of the other. It makes tumor3 homozygous at every tumour suppressor at once, and explains how the one tumor carrying **no** engineered edits still became the most aneuploid (§6.6).
6. **Study A tumors are aneuploid, each with a distinct karyotype.** The three tumors carry different chromosome-scale gains/losses, and no aneuploidy is shared in the same direction by any two — independent clonal evolution. The parental and edited-but-pre-tumor cells remain diploid.
7. **Study B tissues are genomically stable — no aneuploidy.** All six Lats1/2-flox tissues (both genotypes, all ages) show a flat, diploid copy-number profile genome-wide (autosomal copy number 1.94–2.08). The reported oviduct phenotype is **not** driven by large-scale chromosomal instability or aneuploidy. **However, this excludes one mechanism, not all of them:** the two hypotheses ranked most likely at the design stage — iHPV insertional mutagenesis (integration site unresolved) and leaky loxP recombination (assay not yet valid) — **remain untested**, so this is not a "no genomic cause" result. Status is set out in §6.10.
8. **The phenotype affects both strains, but only one carries the transgene** — so an iHPV integration cannot be the whole explanation. What the two strains share is the floxed Lats1/2 alleles and a non-C57BL/6J background; both are live candidates, and the first (a hypomorphic floxed allele) is invisible to WGS and needs RNA or protein to test (§6.10).
9. **The iHPV transgene is present specifically in L1L2H mice.** HPV16 E6/E7 and luciferase reads are detected in all three L1L2H samples and are entirely absent from all three L1L2 samples (perfect specificity), at levels consistent with a fixed germline transgene.
10. **Study B is not congenically pure C57BL/6J.** Each tissue carries ~5–6 million variants against GRCm39 — a genome-wide non-6J inbred background, most likely 129-derived from the ES-cell engineering. This is inherited background, not accumulated mutation, and it shapes how de-novo candidates must be mined (§6.9).

---

## 3. Sample Information

12 samples, paired-end 150 bp WGS (NovaSeq X Plus), aligned to GRCm39 (GENCODE vM35) with nf-core/sarek 3.8.1 (bwa-mem2). Mean autosomal depth in parentheses.

| #  | Label     | Study |  Type  | Group  | Genotype / role                                      |
| :- | :-------- | :---: | :----: | :----: | :--------------------------------------------------- |
| 1  | RO_origin |   A   |  Cell  | parent | Trp53⁺/⁻; Cas9 transgenic; unedited — **matched normal** (~20×) |
| 2  | RO_B1TP   |   A   |  Cell  | edited | Parent cells electroporated with **Brca1 + Pten** sgRNAs (~24×) |
| 3  | RO_B2TP   |   A   |  Cell  | edited | Parent cells electroporated with **Brca2 + Pten** sgRNAs (~24×) |
| 4  | RO_tumor1 |   A   |  Cell  | tumor  | Cell line from a solid tumor grown after injecting B1TP/B2TP cells into a mouse (~26×) — **resolved here to the B1TP lineage** (§6.3) |
| 5  | RO_tumor2 |   A   |  Cell  | tumor  | Cell line from a **different** such tumor (~19×) — **resolved here to the B1TP lineage** (§6.3) |
| 6  | RO_tumor3 |   A   |  Cell  | tumor  | Cell line from a **third** such tumor (~30×) — carries **none** of the intended edits; **resolved here to an un-edited subclone of the parent** (§6.3–6.4) |
| 7  | L1L2_3M   |   B   | Tissue |  L1L2  | Lats1/2 flox, 3 months (~21×)                        |
| 8  | L1L2H_3M  |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 3 months (~20×)                 |
| 9  | L1L2_12M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 12 months (~32×)                       |
| 10 | L1L2H_12M |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 12 months (~28×)                |
| 11 | L1L2_18M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 18 months (~25×)                       |
| 12 | L1L2H_18M |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 18 months (~21×)                |

Total input: 12 samples, gzip FASTQ 462 GiB (measured on disk 2026-07-12).

**Study A design.** All six Study A samples descend from one parental primary-cell line (RO_origin) taken from a **Trp53⁺/⁻; Cas9** transgenic mouse — it already carries one Trp53-null allele and expresses Cas9 constitutively, and it was **not** edited in vitro. Because Cas9 is already present in the cells, editing required delivering **sgRNA only, by electroporation** — no viral or plasmid vector is involved in Study A, and no donor template was used (knockout by non-homologous end joining at the cut sites). Since RO_origin already carries the Trp53 and Cas9 features, using it as the matched normal subtracts them correctly, leaving only the sgRNA-induced edits and tumor-acquired changes. **Tumors 1–3 are three separate solid tumors** grown after injecting the edited cells into mice, each digested into its own cell line; the client did not record which injected line each tumor came from, and resolving that from the WGS is objective A2 (answered in §6.3). *(The iHPV construct belongs to Study B only — it plays no part in Study A.)*

### 3.1 CRISPR guides (CRISPRevolution sgRNA EZ Kit, 1.5 nmol each, Modified)

Nine guides were supplied as images and converted RNA→DNA. **Every guide was independently verified before use:** each maps to a **single unique site** inside its intended target gene, on the strand its name specifies where one is given, and each carries a **perfect NGG PAM** — confirming the guides are genuine and correctly transcribed. Cut site = 3 bp upstream of the PAM.

| Target | Guide ID (as supplied) | Spacer (DNA, 5'→3')    | GRCm39 location                | Strand | PAM | Predicted cut     | Used in     |
| :----- | :--------------------- | :--------------------- | :----------------------------: | :----: | :-: | :---------------: | :---------: |
| Pten   | Pten-32799878          | `GGTGGGTTATGGTCTTCAAA` | chr19:32,777,275–32,777,294    |   −    | AGG | chr19:32,777,278  | B1TP + B2TP |
| Pten   | Pten-32799895          | `TGATAAGTTCTAGCTGTGGT` | chr19:32,777,292–32,777,311    |   −    | GGG | chr19:32,777,295  | B1TP + B2TP |
| Pten   | Pten-32799899          | `GGTTTGATAAGTTCTAGCTG` | chr19:32,777,296–32,777,315    |   −    | TGG | chr19:32,777,299  | B1TP + B2TP |
| Brca1  | (guide 1)              | `GGTTCCGGTAGCCCACGCTC` | chr11:101,422,890–101,422,909  |   +    | TGG | chr11:101,422,906 | B1TP        |
| Brca1  | (guide 2)              | `GGCGTCGATCATCCAGAGCG` | chr11:101,422,905–101,422,924  |   −    | TGG | chr11:101,422,908 | B1TP        |
| Brca1  | (guide 3)              | `TTCTTGTGAGCGTTTGAATG` | chr11:101,422,929–101,422,948  |   −    | AGG | chr11:101,422,932 | B1TP        |
| Brca2  | Brca2+150529497        | `GATAAGCCTCAATTGGTTTG` | chr5:150,452,945–150,452,964   |   +    | AGG | chr5:150,452,961  | B2TP        |
| Brca2  | Brca2−150529492        | `AAAGCTCCTCAAACCAATTG` | chr5:150,452,954–150,452,973   |   −    | AGG | chr5:150,452,957  | B2TP        |
| Brca2  | Brca2−150529524        | `AGGTTCAGAATTGTATGGGG` | chr5:150,452,986–150,453,005   |   −    | GGG | chr5:150,452,989  | B2TP        |

**All nine cut sites fall inside a coding exon of the Ensembl-canonical transcript of their target** — Pten CDS exon 5 (`ENSMUST00000249247.1`), Brca1 CDS exon 6 (`ENSMUST00000017290.11`), Brca2 CDS exon 3 (`ENSMUST00000044620.11`, *Brca2-201*, CCDS39411.1) — so a disruptive indel at any of them is expected to be loss-of-function rather than silent.

**All three targets use the same multi-guide strategy:** the three guides per gene are clustered within a narrow window (Pten 32,777,278–32,777,299 = 22 bp; Brca1 101,422,906–101,422,932 = 27 bp; Brca2 150,452,957–150,452,989 = 33 bp), with overlapping guides cutting the same point from opposite strands. This predicts either small indels at a cut site or excision of the fragment between the outermost cuts — the latter is exactly what is observed in B2TP (§6.3).

**Note on the guide IDs:** the coordinates embedded in the product names **are not GRCm39 positions**. `Brca2+150529497` sits ~76.5 kb from where that guide actually maps on GRCm39, and `Pten-32799878` ~22.6 kb from its true position — so using these numbers directly against GRCm39 retrieves the wrong locus. That the offset **differs per locus** (76.5 kb vs 22.6 kb) rules out a simple constant shift and is what an assembly coordinate change produces; the numbering is most consistent with GRCm38/mm10, though this was not verified directly (no GRCm38 reference was used here). Every guide was located by sequence, so no result is affected.

---

## 4. Analysis Rationale and Decision Criteria

| Question                        | Approach & why it works without a dedicated wild-type control                                                                       | Threshold / criterion                                                      |
| :------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------- |
| Copy number / aneuploidy        | Ratio of genome-internal coverage → copy number; reference-free. GRCm39 = C57BL/6J ≈ the animals' own normal.                     | Chromosome-median CN ≥2.5 = gain, ≤1.5 = loss                            |
| Edit verification (A1)          | Read directly at each sgRNA cut site; compare edited/tumor vs the RO_origin parent. **Read/CIGAR level, not caller output** — see §6.3 | Wild-type reads absent ⇒ biallelic knockout; indel present in edited sample, absent (0/0) in RO_origin |
| Soft-clipped reads as evidence  | This locus is intrinsically clip-prone (~30% of reads clip even in the unedited parent), so clip *fraction* is not specific         | A clip counts only if ≥20 bp **and** its breakpoint **stacks recurrently on a cut-site base** |
| Negative-control gate           | The unedited parent must score as unedited, or the method is wrong                                                                  | RO_origin must show ~0 editing (asserted in code)                          |
| Tumor lineage (A2)              | Brca1-KO ⇒ B1TP lineage; Brca2-KO ⇒ B2TP lineage. A line with **no** wild-type allele cannot give rise to a wild-type descendant | Presence/absence of cut-site indel; wild-type read fraction               |
| Sample identity (A2)            | Naive genotype concordance is confounded by tumor LOH (het→hom), which would penalise the most aneuploid sample; **private alleles** are immune, since LOH removes alleles but cannot create one the parent never had | Homozygous private alleles ≈ 0 ⇒ same animal; thousands ⇒ different animal |
| iHPV presence / integration (B) | Transgene reads are unmapped against plain GRCm39 → extract and align to HPV16/EGFP/Luc markers; L1L2 = internal negative control. | Construct-mapping reads present in L1L2H, absent in L1L2                   |
| De-novo candidates (B)          | Subtract known mouse strain background (Sanger Mouse Genomes Project [MGP], all strains) from each sample's germline calls.          | Variant private to sample (not in MGP), recurrent or genotype-differential |

---

## 5. Methods

- **Alignment / preprocessing:** nf-core/sarek 3.8.1, bwa-mem2, GRCm39 (GENCODE vM35), duplicate marking, BQSR skipped.
- **Variant calling:** Study A somatic = GATK Mutect2 (each tumor vs RO_origin) + TIDDIT SV; Study B germline = GATK HaplotypeCaller + CNNScoreVariants + TIDDIT SV.
- **Copy number:** mosdepth 500 kb binned depth → normalised to per-sample autosomal median → copy number = 2 × ratio; chromosome-level medians and cohort heatmap (mouse chromosomes are acrocentric, so whole-chromosome CN ≈ arm-level).
- **Guide localisation:** each spacer located on GRCm39 by sequence, both strands (`seqkit locate`); PAM confirmed from the 3 bp immediately 3' of each protospacer; cut site = 3 bp upstream of the PAM.
- **Edit verification:** multi-sample `bcftools mpileup`/`call` at cut-site windows (MAPQ≥20, BAQ≥15), **with read-level allele quantification as the primary read-out**. Standard variant callers under-report at CRISPR cut sites, where edited reads carry large deletions or are heavily soft-clipped rather than presenting as simple indels; genotypes here are therefore taken from the reads themselves. Read-level analysis (custom pysam): reads with MAPQ≥20, duplicates excluded, classified by CIGAR into wild-type (spanning the cut-site core with ≥25 bp flanking margin, no indel/clip), the excision allele, other cut-site indels, or cut-site soft-clips (≥20 bp, breakpoint inside the core); clip breakpoints tabulated per base to separate recurrent from background clipping.
- **Sample identity fingerprint:** joint `bcftools mpileup/call` across all six Study A samples over 38 windows × 500 kb (19.0 Mb) spanning all 19 autosomes; biallelic SNPs retained. A "private allele" = sample carries ALT with ≥3 supporting reads while RO_origin has **0** ALT reads at ≥10× depth. Results broken down per chromosome and per allele fraction.
- **iHPV integration:** unmapped reads extracted per sample and aligned (bwa-mem2) to a marker reference (HPV16 NC_001526.4 + EGFP + firefly luciferase); genomic anchors from mapped reads with unmapped mates.
- **De-novo candidates:** Sanger Mouse Genomes Project v8 (REL-2021, GRCm39) SNP+indel VCFs subtracted from each sample's germline calls with `bcftools isec` (chromosome names reconciled chr1↔1).

---

## 6. Results

### 6.1 Sequencing data quality — all 12 samples pass

Every sample is of good and comparable quality; no sample was excluded, and no result below is limited by data quality.

| Sample    | Mapped % | Duplicate % | Mean depth | Genome ≥10× | Insert size | Mean base quality | Error rate |
| :-------- | :------: | :---------: | :--------: | :---------: | :---------: | :---------------: | :--------: |
| RO_origin |  99.99   |    10.78    |   19.4×    |     86%     |    139 bp   |       38.7        |  4.5×10⁻³  |
| RO_B1TP   |  99.99   |    14.22    |   23.2×    |     89%     |    123 bp   |       38.8        |  4.6×10⁻³  |
| RO_B2TP   |  99.99   |    14.29    |   23.6×    |     89%     |    128 bp   |       38.7        |  4.8×10⁻³  |
| RO_tumor1 |  99.99   |    11.83    |   26.0×    |     90%     |    136 bp   |       38.7        |  4.6×10⁻³  |
| RO_tumor2 |  99.99   |    10.96    |   19.2×    |     84%     |    131 bp   |       38.7        |  4.6×10⁻³  |
| RO_tumor3 |  99.99   |    14.08    |   31.6×    |     91%     |    125 bp   |       38.4        |  4.8×10⁻³  |
| L1L2_3M   |  99.99   |     9.43    |   21.5×    |     88%     |    141 bp   |       38.8        |  4.7×10⁻³  |
| L1L2H_3M  |  99.99   |    11.83    |   20.1×    |     86%     |    132 bp   |       38.5        |  4.7×10⁻³  |
| L1L2_12M  |  99.99   |     7.68    |   31.9×    |     92%     |    140 bp   |       38.6        |  4.7×10⁻³  |
| L1L2H_12M |  99.99   |    10.70    |   27.6×    |     92%     |    137 bp   |       38.9        |  4.7×10⁻³  |
| L1L2_18M  |  99.99   |     9.30    |   24.8×    |     90%     |    139 bp   |       38.7        |  4.8×10⁻³  |
| L1L2H_18M |  99.99   |    10.30    |   20.5×    |     87%     |    137 bp   |       38.7        |  5.1×10⁻³  |

- **Mapping is essentially complete** — 99.99% of QC-passed reads align to GRCm39 in every sample (e.g. RO_origin: 53,411 unmapped of 458.9 M). Upstream, fastp retained 99.0% of raw reads, so ~99% of all sequenced reads are used.
- **Duplicate rates are low** (7.7–14.3%) and **base quality is high** (mean Q38.4–38.9; Q30 rate 93.8%), with a uniform error rate of ~0.5%.
- **Depth is 19.2–31.9×**, with 84–92% of the genome covered at ≥10×. This is ample for the copy-number, edit-verification and transgene analyses reported here. It is at the lower end for *sensitive* somatic point-mutation calling, which is one reason the somatic burden is not quoted (§6.11).
- **One library characteristic worth noting:** the mean insert size (123–141 bp) is **shorter than the 150 bp read length** in all 12 samples, so read pairs overlap substantially. This is systematic across the batch rather than a per-sample defect. It has two consequences: the effective independent coverage is somewhat below the nominal depth (overlapping mates re-sequence the same molecule), and paired-end SV detection is less sensitive for a given depth — relevant to the SV counts in §6.8, and a reason those are treated as provisional.

Full per-sample reports: `qc/multiqc_studyA.html`, `qc/multiqc_studyB.html`.

### 6.2 Copy number & aneuploidy (headline)

Chromosome-level copy number (autosomes) across the cohort:

- **Study B (all six Lats1/2 tissues):** flat and diploid — autosomal copy number ranges only **1.94–2.08** in every sample, at every age, in both genotypes. No chromosome-scale gain or loss.
- **Study A parent + edited cells (RO_origin, B1TP, B2TP):** also flat/diploid (1.91–2.14).
- **Study A tumors:** clearly aneuploid, each distinct (chromosome-median copy number in brackets; calls at CN ≥2.5 = gain, ≤1.5 = loss — full table in `cnv_ploidy/aneuploidy_calls.tsv`) —
  - **tumor1:** gains chr8 [2.60], chr10 [2.69], chr11 [2.52]; loss chr18 [1.46]
  - **tumor2:** gain chr5 [2.52] (plus milder shifts)
  - **tumor3:** losses chr8 [1.43], chr12 [1.36], chr13 [1.43], chr14 [1.37]; gain chr15 [2.51] — **five affected chromosomes, the most-rearranged genome**

**No aneuploidy is shared in the same direction by any two tumors** — the karyotypes do not overlap. The one chromosome altered in two tumors is **chr8, and in opposite directions** (gained in tumor1 [2.60], lost in tumor3 [1.43]), which reinforces rather than weakens the point. These are three independent clonal events, not one tumor sampled three times — consistent with their distinct lineages (§6.3: tumors 1 and 2 from B1TP, tumor3 from an un-edited subclone).

**Interpretation.** The Study A tumors act as an internal positive control: the pipeline detects clear aneuploidy where it exists. Against that, the uniformly flat Study B profiles are a genuine biological result — **the Lats1/2-flox tissues have stable, diploid genomes with no aneuploidy, so the oviduct phenotype is not explained by large-scale genomic instability.** This is consistent with the engineered elements being un-activated (Lats1/2 not deleted). See `cnv_ploidy/cohort_cn_heatmap.png` and per-sample `*.cn_profile.png`.

### 6.3 CRISPR edit verification (Study A)

All nine cut sites were genotyped at read level (§5). Cut-site genotypes across the three target genes — RO_origin is homozygous reference throughout, confirming these are edits, not background:

| Gene / site             | RO_origin |     B1TP      |          B2TP          |         tumor1         |     tumor2      | tumor3 |
| :---------------------- | :-------: | :-----------: | :--------------------: | :--------------------: | :-------------: | :----: |
| Brca1 (chr11:101422906) |    0/0    | **indel**     |          0/0           | **indel (biallelic)**  | **indel (hom)** |  0/0   |
| Brca2 (chr5:150452957–989) |  wild-type |   wild-type   | **biallelic KO**  |       wild-type        |    wild-type    | wild-type |
| Pten (chr19:32777294)   |    0/0    | **indel**     | **indel (hom)**        |    **indel (hom)**     |      indel      |  0/0   |

**Brca2 allele composition at the cut site.** Counts are **independent DNA fragments**, not reads: this library's inserts (123–141 bp) are shorter than the 150 bp reads (§6.1), so the two mates of a pair cover the same molecule and counting reads would double-count the evidence.

| Sample        | Informative fragments | **Wild-type** | 31 bp excision | Disrupted at cut site | Other |
| :------------ | :-------------------: | :-----------: | :------------: | :-------------------: | :---- |
| RO_origin     |           5           |       4       |       0        |          1            | –     |
| RO_B1TP       |           6           |       5       |       0        |          1            | –     |
| **RO_B2TP**   |         **9**         |     **0**     |       1        |        **7**          | 1 fragment with a 3 bp deletion |
| RO_tumor1     |          15           |      11       |       0        |          4            | –     |
| RO_tumor2     |           9           |       5       |       0        |          4            | –     |
| RO_tumor3     |          12           |       4       |       0        |          7            | 1 fragment with a 1 bp deletion |

- **B1TP:** Brca1 KO + Pten KO — matches its Brca1+Pten design.
- **B2TP — Brca2 knockout confirmed, biallelic.** Two independent lines of evidence establish this:
  - **No wild-type Brca2 allele survives.** B2TP is the only sample with **0 wild-type fragments** (0 of 9), while every other sample retains a clear wild-type population (4–11 fragments). Were one Brca2 allele still wild-type, roughly half of those 9 fragments should have been wild-type; observing none has a probability of about 0.002.
  - **Coverage collapses at the cut site and nowhere else.** B2TP's Brca2 coverage is normal across the gene (23.2×) but falls to **0.57 of its own baseline** precisely at the cut window (13.3×), the expected consequence of edited molecules failing to align cleanly. Every other sample sits at 0.87–1.23. This measure does not depend on read counting at all.
- **How each allele is disrupted is not resolved, and the data do not support naming them.** Effective depth at the cut window is only ~7×, giving 9 informative fragments in total — enough to establish that no wild-type allele remains, but not to reconstruct two allele structures. Specifically:
  - **The predicted excision is observed, on a single molecule.** One fragment carries a 31 bp deletion at chr5:150,452,958–150,452,988, removing `TGGTTTGAGGAGCTTTCCTCAGAAGCCCCCC`. Its endpoints are not merely near the guides — the 5' end abuts the overlapping guide pair's cut sites (150,452,957/961) and the 3' end abuts the third guide's cut (150,452,989), i.e. exactly the fragment predicted to be excised when all three guides cut. Falling wholly within CDS exon 3 from CDS nucleotide 91 (**codon ~31 of 3,329**), and not being a multiple of 3, it would cause a frameshift and premature termination in the first 1% of the protein. **On one molecule this is a consistent observation, not an established allele** — it should be read as the excision having occurred in this cell line, not as a genotype.
  - **Most disrupted molecules break at the third guide's cut.** Seven of the nine fragments are soft-clipped, and their breakpoints **stack on a single base** (150,452,989/990) rather than scattering as in every other sample. The clipped sequence maps only ambiguously elsewhere (MAPQ 0, repetitive), so whether this represents a complex indel, an insertion, or a larger rearrangement is not resolvable with short reads.
  - A single fragment carries a 3 bp deletion at the third cut; one molecule is not evidence of an allele and it is not interpreted here.
- **B2TP = Brca2 + Pten knockout** is therefore established by direct observation that Brca2 is disrupted on both alleles, rather than inferred from its Pten/Brca1 pattern. Resolving the exact lesion on each allele would require targeted deep sequencing (amplicon or long-read) across the cut site — not more WGS.
- **tumor1, tumor2:** carry both Brca1 and Pten indels → **B1TP (Brca1+Pten) lineage**; both retain wild-type Brca2, as expected.
- **tumor3:** wild-type at **all three** targeted genes. Its only Brca2 indel call — a 1 bp deletion at 150,452,983 on a single fragment — falls inside a **7-C homopolymer** (`AAG`**`CCCCCCC`**`ATACAATTCTG`, 150,452,983–150,452,989), the classic context for a sequencing artifact, and is discounted. The lineage logic is decisive: **B2TP has no wild-type Brca2 allele**, so no descendant of B2TP could regain one; tumor3 is ~47% wild-type at this locus and therefore **cannot descend from B2TP**. Pten is equally decisive in the other direction — **both** injected lines carry a Pten edit (B2TP homozygously), so descent from either line would leave a detectable Pten lesion, and tumor3 has none (Pten sites 17–18× reads, unambiguous reference — a true absence of edit, not a coverage gap). Tumor3 descends from neither edited line.

### 6.4 Tumor3 identity — same animal as the parent

Tumor3 being wild-type at all three targeted genes leaves a question of **identity**: an un-edited escaper subclone of the parent, or a mis-tracked sample from a different mouse? B1TP/B2TP/tumor1/tumor2 — all known to descend from origin — provide the "same animal" baseline. 19,032 evaluable sites.

| Sample        | Private alleles (raw) | Excluding chr3/chr19 | **Homozygous private (excl.)** |
| :------------ | :-------------------: | :------------------: | :----------------------------: |
| RO_B1TP       |          245          |         223          |               0                |
| RO_B2TP       |          227          |         216          |               0                |
| RO_tumor1     |          264          |         231          |               0                |
| RO_tumor2     |          219          |         182          |               1                |
| **RO_tumor3** |      **2,851**        |         644          |            **8**               |

Tumor3's raw count (2,851, ~11× the other samples) does not indicate a different animal. It resolves into two components, neither of them inherited variation:

- **By chromosome:** 1,154 of tumor3's 1,218 homozygous private alleles fall in a **single 500 kb window on chr3 (~47.9–48.4 Mb)** and 56 on chr19; **all 17 remaining autosomes contribute zero**. A different animal differs on *every* chromosome; this is a local event. That chr3 window also shows anomalous coverage in tumor3 alone (0.35× its own baseline, vs 0.69–1.00 in every other sample) — a tumor-specific structural loss whose residual mismapped reads generate spurious homozygous calls.
- **By allele fraction:** excluding those two regions, tumor3 retains just **8 homozygous private alleles**. Its remaining excess sits at **VAF 0.3–0.7** — the signature of **clonal somatic mutation**, expected in the most rearranged of the three tumors. Inherited variation would instead pile up at VAF ~1.0.

**Conclusion:** tumor3 shares RO_origin's germline genome. It is the same animal, and a swap with a different mouse is excluded. With §6.3, **tumor3 arose from an un-edited (editing-escaped) subclone of the Trp53⁺/⁻;Cas9 parent** — and its heavy aneuploidy therefore developed independently of any Brca/Pten editing.

**One caveat, stated plainly:** this method distinguishes *different genomes*, not *different individuals of an inbred line*. Littermates of an inbred colony are near-genetically identical, so a mix-up between two such animals would be invisible to any fingerprint. What is established is that tumor3's genome is the parent's genome — precisely what the un-edited-subclone explanation requires, and what a swap with an unrelated or differently-engineered animal would have broken.

### 6.5 Loss of heterozygosity — all three tumors lost Trp53, by three different mechanisms

The parent is **Trp53⁺/⁻** by design, so the classic route to a Trp53-null tumor is a "second hit" that removes the remaining functional allele. Coverage cannot see this: Trp53 read depth is flat and unchanged in all six samples, including the parent, meaning the engineered null allele is a small lesion rather than a deletion. **Loss of heterozygosity answers it instead, and requires no knowledge of the allele's design** — we simply ask whether SNPs that are heterozygous in RO_origin have become homozygous in each descendant. B1TP and B2TP — edited but never passaged through a tumor — provide the negative control and define the noise floor.

**Genome-wide LOH (11,272 origin-heterozygous SNPs across all 19 autosomes):**

| Sample | Genome-wide LOH |
| :--- | :---: |
| RO_B1TP | 5% |
| RO_B2TP | 4% |
| RO_tumor1 | 10% |
| RO_tumor2 | 13% |
| **RO_tumor3** | **95%** |

**chr11 scan across the Trp53 locus (Trp53 = 69.47 Mb):**

| chr11 window | Informative SNPs | B1TP | B2TP | **tumor1** | tumor2 | tumor3 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| 55–56 Mb | 367 | 3% | 3% | **3%** | 94% | 95% |
| 65–66 Mb | 634 | 1% | 1% | **1%** | 98% | 97% |
| **69–70 Mb (Trp53)** | 725 | 2% | 2% | **92%** | 93% | 94% |
| 72–73 Mb | 680 | 2% | 2% | **3%** | 96% | 97% |

**All three tumors are homozygous at Trp53; none of the pre-tumor cells are.** The mechanism differs in each, which is itself the point — three independent routes converging on the same locus:

- **tumor1 — a focal LOH tract centred on Trp53.** This is the decisive observation. Tumor1's genome-wide LOH is only 10%, and on chr11 it is heterozygous **1% at 65–66 Mb and 3% at 72–73 Mb** — yet **92% at 69–70 Mb, where Trp53 sits**. Heterozygosity is intact 3 Mb away on *both* sides. A random event does not carve a narrow homozygous tract precisely over one tumour-suppressor gene and stop; this is selection for Trp53 loss, written directly into the data.
- **tumor2 — a large chr11 regional LOH** (94–98% from 55 Mb through 73 Mb) that includes Trp53, against a 13% genome-wide baseline.
- **tumor3 — genome-wide LOH** (§6.8); its Trp53 homozygosity comes with the whole genome rather than being locally targeted.

**What this does and does not establish.** Each tumor retains only one Trp53 haplotype, so each is homozygous — either null/null or wild-type/wild-type. Distinguishing which haplotype survived would require the design of the engineered "−" allele (what the lesion actually is), which we do not have; with it, the existing data are sufficient to genotype it directly. But the inference is not finely balanced: retaining the wild-type copy would *restore* p53 function and confer no advantage, whereas losing it is the canonical second hit in Trp53⁺/⁻ tumors. Three independent tumors converging on Trp53 homozygosity — one of them focally — is what selection for **Trp53 loss** looks like, and it supplies the driver that the Brca/Pten edits do not explain in tumor3 (§7.1).

### 6.6 Tumor3 has undergone genome-wide, copy-neutral LOH

Tumor3's LOH is not confined to Trp53 or to any chromosome: **95% of origin-heterozygous sites genome-wide are homozygous in tumor3**, versus 4–5% in the edited pre-tumor cells. The allele-fraction spectrum settles it beyond any genotype-calling artefact — at sites where the parent is a clean heterozygote, tumor3's reads are **bimodal at 0 and 1** with almost nothing at 0.5:

| ALT allele fraction | RO_origin | RO_B1TP | RO_B2TP | RO_tumor1 | RO_tumor2 | **RO_tumor3** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.0–0.1 | 2.0% | 5.1% | 4.6% | 9.2% | 9.6% | **40.7%** |
| 0.4–0.6 (heterozygous) | 50.7% | 49.1% | 52.2% | 19.3% | 38.1% | **1.0%** |
| 0.9–1.0 | 0.0% | 0.2% | 0.1% | 1.6% | 2.6% | **53.2%** |
| Sites assessed | 3,433 | 1,711 | 1,865 | 2,430 | 1,010 | 2,528 |

Tumor3 has essentially no heterozygous sites left. Crucially, **its copy number is still diploid across most of the genome** — chr2 2.06, chr7 2.03, chr16 2.09, chr18 2.19 — while those same chromosomes are 94–99% homozygous. Homozygosity at normal copy number is **copy-neutral LOH**, and genome-wide copy-neutral LOH has one standard explanation: **loss of one entire haploid genome followed by duplication of the remainder** (haploidisation and endoreduplication). The aneuploidies in §6.2 (chr8, chr12, chr13, chr14 losses; chr15 gain) are superimposed on that event.

This is the most extreme genomic phenotype in the cohort, and it belongs to **the one tumor that carries none of the intended edits** (§6.3). It hands tumor3 homozygosity at *every* tumour suppressor at once, Trp53 included — a single catastrophic route to the same endpoint the Brca/Pten edits were designed to engineer.

It also sets the limit on how tumor3's identity can be assessed: at 95% homozygosity, any test based on genotype concordance with the parent would score tumor3 as broadly discordant regardless of its true origin. The private-allele metric in §6.4 is unaffected, since LOH removes alleles but cannot create ones the parent never carried.

### 6.7 iHPV transgene detection (Study B)

Construct-marker reads (unmapped reads aligned to HPV16 + EGFP + luciferase):

| Sample    | Construct reads | HPV16 | EGFP | Luciferase |
| :-------- | :-------------: | :---: | :--: | :--------: |
| L1L2_3M   |        0        |   0   |  0   |     0      |
| L1L2H_3M  |       45        |  16   |  0*  |     71     |
| L1L2_12M  |        0        |   0   |  0   |     0      |
| L1L2H_12M |       83        |  41   |  0*  |    123     |
| L1L2_18M  |        0        |   0   |  0   |     0      |
| L1L2H_18M |       58        |  37   |  0*  |     75     |

**The iHPV transgene is present in all three L1L2H samples and absent in all three L1L2 samples** — a perfectly specific result. Read counts do not increase with age, consistent with a fixed germline transgene rather than an age-accumulating event.

**\* The zero EGFP counts carry no biological information and must not be read as the EGFP cassette being absent.** The GFP sequence used in the marker reference encodes the correct protein but is a low-GC codon variant (714 bp, 40.3% GC), whereas a mammalian CAG-driven construct carries human-codon-optimised EGFP (~720 bp, 61–62% GC). The two differ at roughly 30% of nucleotide positions — beyond what short-read alignment tolerates — so EGFP-derived reads cannot map to this reference even when the cassette is physically present. The consequence matters: the EGFP cassette is the **stop** element of the lox-stop-lox, and its integrity is precisely the test for whether E6/E7 has been de-repressed (§6.10). **That test remains outstanding**; this table is uninformative about it in either direction. It requires the actual EGFP sequence from the vector map.

### 6.8 Structural variants

TIDDIT SV calls, PASS-filtered (`sv/sv_counts.tsv`):

| Sample    | Study | Total | PASS | Deletions | Inversions | Breakends |
| :-------- | :---: | ----: | ---: | --------: | ---------: | --------: |
| RO_origin |   A   |  6,073 | 3,803 |  2,869 |  25 |   540 |
| RO_B1TP   |   A   |  5,417 | 3,725 |  3,001 |  15 |   400 |
| RO_B2TP   |   A   |  5,925 | 3,857 |  3,114 |  14 |   358 |
| RO_tumor1 |   A   |  6,444 | 3,880 |  2,993 |  17 |   426 |
| RO_tumor2 |   A   |  4,724 | 3,362 |  2,694 |  16 |   368 |
| **RO_tumor3** | A | **10,118** | **5,416** | **4,347** | **27** | 508 |
| L1L2_3M   |   B   | 16,393 | 10,425 |  8,261 |  71 | 1,250 |
| L1L2H_3M  |   B   | 13,201 |  8,971 |  7,181 |  75 | 1,040 |
| L1L2_12M  |   B   | 27,337 | 14,516 | 11,860 |  77 | 1,168 |
| L1L2H_12M |   B   | 19,299 | 11,331 |  9,103 |  67 | 1,162 |
| L1L2_18M  |   B   | 21,038 | 12,769 | 10,269 |  75 | 1,364 |
| L1L2H_18M |   B   | 16,581 | 11,041 |  8,926 |  65 | 1,186 |

**These are raw per-sample calls, not somatic SVs, and the absolute numbers should not be read as SV burden.** They are counts against the reference and therefore include every inherited SV a sample carries. The Study B counts (9,000–14,500 PASS) are two-to-three times the Study A counts for exactly this reason — those tissues sit on a non-6J strain background (§6.9), so most of their "SVs" are inherited strain differences from GRCm39, not lesions. Comparisons are only meaningful *within* an arm, against a matched background.

**Within Study A, one signal survives that caveat: tumor3 is the clear outlier** — 5,416 PASS calls against 3,362–3,880 for the parent, both edited lines and the other two tumors (~1.4× the highest, ~1.6× the lowest), driven mainly by deletions (4,347 vs 2,694–3,114). All six Study A samples share the same genetic background, so this excess is not a background artefact. It is independent corroboration of the copy-number result (§6.2): **tumor3 has the most rearranged genome of the three tumors** — and it is the tumor carrying none of the intended edits (§6.3).

Two limits on interpretation: TIDDIT counts here are per-sample rather than tumor-minus-normal (somatic SV calling against RO_origin is a listed next step, §7.3), and the short insert size (§6.1) reduces paired-end SV sensitivity uniformly across samples. Neither affects the relative comparison within Study A.

### 6.9 De-novo candidate variants (Study B)

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

- Each sample carries **~5–6 million variants relative to GRCm39** — a firm and informative result. Because **GRCm39 is the C57BL/6J ("6J") reference assembly**, a genuinely pure C57BL/6J animal sits only tens of thousands of variants from it (ordinary colony drift). A count in the *millions* is therefore a clear signal, not noise: **these tissues carry a genome-wide non-6J classical-inbred-strain background — i.e. the line is not congenically pure C57BL/6J.** This is inherited strain background, **not** mutation accumulated since the line was made: breeding over many generations, time in culture, and the gene-editing itself together add at most tens of thousands of *new, sample-private* variants (two to three orders of magnitude too few for 5–6 M), and the empirical ceiling for same-strain drift — C57BL/6J vs C57BL/6N, separated for ~70 years — is likewise only tens of thousands. Decisively, **~90% of the load here (≈4.5–5.4 M of the 5–6 M) is catalogued in known MGP strains** and removed by subtraction; a variant can be in MGP only if it is an inherited strain polymorphism, whereas culture/breeding/editing mutations would be sample-private and would *not* subtract out.
- **Most likely source (educated guess): a 129-derived engineering background.** Gene-targeted and floxed alleles are typically generated in **129-derived embryonic stem cells**, and a 129 background is carried genome-wide unless the line is fully backcrossed to congenic purity on C57BL/6J. This line's ES-cell provenance is not documented in the materials provided, so we present 129 as the leading explanation rather than a certainty — but it is well supported: 129 is both the standard ES-cell donor for this kind of allele and a quantitative match, since **the 129 strain lies ~4–6 million SNPs/indels from the 6J reference — essentially the *whole* of the ~5–6 M observed here.** That near-complete overlap means the background is not a small linked remnant around the engineered loci but a **substantial, genome-wide 129-type genome**. A definitive call on the exact donor strain — and on whether it is pure 129 or a 129×B6 mix — is available on request via a dedicated strain-assignment analysis.
- After MGP subtraction, **~0.57–0.66 million private variants** remain per sample, distributed uniformly across all chromosomes (proportional to chromosome length) — i.e. residual background divergence not covered by MGP, not localised de-novo events.

Genuine causal de-novo events number in the dozens-to-hundreds, not hundreds of thousands. **Identifying them from this background requires functional restriction to high-impact coding consequences (frameshift / stop-gain / splice, via VEP or snpEff) combined with recurrence and L1L2-vs-L1L2H differential filtering.** That functional-annotation step is the recommended immediate next analysis; the per-sample private call sets are provided for it (`candidates_denovo/*.private.vcf.gz`). Critically, this does not affect the primary Study B conclusion: the copy-number analysis (§6.2) already establishes these genomes are structurally stable and diploid.

### 6.10 The oviduct phenotype — what has and has not been tested

This is the question Study B exists to answer, so we state its status plainly rather than leaving it to be inferred from the sections above. Three mechanisms could produce an age-dependent oviduct phenotype in mice whose engineered elements are, by design, un-activated. **Only one has been tested.**

| Hypothesis (priority) | Status | What we can say |
| :--- | :--- | :--- |
| **1. iHPV insertional mutagenesis** — the construct integrates into and disrupts a gene | **Not tested** | We established the transgene is present in L1L2H and absent in L1L2 (§6.7). **Where it integrated is not resolved**, so whether it disrupts a gene is unknown. Needs the full vector map to capture junction-spanning reads. |
| **2. Leaky / somatic loxP recombination** — age-dependent Cre-independent excision deleting Lats1/2 exons or de-repressing E6/E7 | **Not tested** | The read-out is whether the EGFP "stop" cassette is intact, and our EGFP marker sequence is the wrong codon variant, so the assay never ran (§6.7). Lats1/2 floxed-exon integrity likewise not assessed. |
| **3. Large-scale genomic instability / aneuploidy** | **Tested — negative** | Flat, diploid copy-number profiles in all six tissues, both genotypes, all ages (§6.2). This mechanism is **excluded**, and that exclusion is solid: the Study A tumors act as a positive control proving the method detects aneuploidy where it exists. |
| **3b. De-novo point/indel variants** | **Incomplete** | Per-sample private variant sets are produced (§6.9) but not yet restricted to high-impact consequences, so no candidate gene list exists. |

**So the honest position is: we have excluded chromosomal instability as the cause, and have not yet tested the two mechanisms that the study design ranked as most likely.** The report's Study B findings should not be read as "no genomic cause was found".

**One structural observation does bear on the biology now, and it constrains hypothesis 1.** You report the abnormality in **both** strains, but only L1L2H carries the iHPV construct — L1L2 has no transgene at all (§6.7, zero construct reads at every age). An iHPV integration therefore **cannot explain the phenotype in L1L2**. Whatever is shared between the two strains is the more promising place to look, and two things are:

- **The floxed Lats1/2 alleles themselves.** A loxP insertion in an intron can reduce expression of the gene it sits in; a floxed allele is not always phenotypically silent before Cre. If either Lats1/2 allele is hypomorphic, both strains carry reduced Hippo-pathway dosage from birth, which is a plausible route to a slow, age-dependent phenotype with no Cre anywhere. **WGS cannot test this** — it needs Lats1/2 RNA or protein from the affected tissue — but it is consistent with everything seen here, including the genomes being stable and diploid.
- **The strain background.** These mice are not congenically pure C57BL/6J; they carry a genome-wide non-6J (most likely 129-type) background (§6.9). If "abnormal" was judged against C57BL/6 expectations, part of the phenotype may be strain characteristic rather than a consequence of the engineering. A wild-type littermate of the same background would separate the two immediately.

### 6.11 Somatic point mutations (Study A) — not yet a usable mutation burden

Mutect2 tumor-vs-origin calling was run for all five pairs (`somatic/somatic_counts.tsv`):

| Pair (vs RO_origin) | PASS total | SNV | Indel |
| :------------------ | ---------: | --: | ----: |
| RO_B1TP             |  90,680 | 49,965 | 22,613 |
| RO_B2TP             |  87,901 | 48,822 | 22,257 |
| RO_tumor1           |  86,549 | 49,886 | 21,772 |
| RO_tumor2           |  79,920 | 47,965 | 18,324 |
| RO_tumor3           | 305,795 | 247,362 | 45,537 |

**These counts are not mutation burdens and should not be quoted as such.** An internal control fixes their scale: **RO_B1TP and RO_B2TP are the parental cells a few passages after sgRNA electroporation**, and biologically differ from RO_origin at essentially nothing beyond their cut sites — a handful of events. They score ~88,000–91,000. That ~90,000-call floor is therefore artefact, arising because this run lacks a contamination estimate and a panel-of-normals, so recurrent sequencing and alignment errors pass as somatic calls.

Against that floor, the only feature that stands out is **tumor3 at ~3.4× above it** (305,795 vs ~80,000–91,000) — directionally consistent with its status as the most rearranged genome (§6.2, §6.8), but not a quantitative result: an artefact floor this large cannot be subtracted reliably, and part of the excess may itself reflect tumor3's aneuploidy and structural anomalies degrading local alignment.

Quoting real burdens requires the filtering listed in §7.3 (VAF/depth thresholds, orientation-bias filtering, a panel-of-normals). Meanwhile the Study A conclusions rest on copy number (§6.2), edit verification (§6.3), identity (§6.4) and SV (§6.8), none of which depend on these counts.

---

## 7. Conclusions

| Conclusion                                                                     | Confidence | Evidence |
| :----------------------------------------------------------------------------- | :--------: | :------- |
| The three Brca2 guides are genuine and target Brca2 exon 3                     |    High    | §3.1 — unique GRCm39 hit each, correct strand, perfect NGG PAM |
| **B2TP carries the intended Brca2 knockout, and it is biallelic**              |    High    | §6.3 — 0 of 9 wild-type fragments (p≈0.002 if one allele were wild-type); plus coverage 0.57× of own baseline confined to the cut window |
| The predicted 31 bp dual-cut excision occurred in B2TP (frameshift, codon ~31/3,329) | Moderate — single molecule | §6.3 — endpoints coincide exactly with the guide cut sites, but only 1 fragment carries it |
| The identity of each individual Brca2 allele                                    | **Not established** | §6.3 — ~7× effective depth / 9 fragments cannot resolve two allele structures; needs amplicon or long-read sequencing |
| Most disrupted B2TP molecules break at the third guide's cut                    | Moderate–High | §6.3 — breakpoints of 7 fragments stack on 150,452,989/990 vs scattered elsewhere |
| **B2TP = Brca2 + Pten knockout**, by direct observation rather than inference   |    High    | §6.3 |
| B1TP = Brca1 + Pten KO confirmed                                                |    High    | §6.3 |
| Tumors 1 & 2 arose from the B1TP lineage; neither carries a Brca2 edit          |    High    | §6.3 |
| **Tumor3 arose from an un-edited subclone of the parent, not B1TP/B2TP**        |    High    | §6.3 (wild-type at all three genes) + §6.4 (shares origin's germline) |
| **Tumor3 is the same animal as RO_origin — sample swap excluded**               | High (see §6.4 caveat re inbred littermates) | §6.4 — 8 homozygous private alleles vs thousands expected for a different mouse |
| Tumor3's chr3 ~47.9–48.4 Mb anomaly is a tumor-specific structural event, not an identity signal | Moderate–High | §6.4 — confined to 1 window; coverage 0.35× of own baseline |
| **All three tumors are homozygous at Trp53 (LOH); the pre-tumor edited cells are not** |    High    | §6.5 — 92–94% LOH at Trp53 vs 2–3% in B1TP/B2TP |
| **Tumor1's Trp53 LOH is focal — direct evidence of selection at Trp53**         |    High    | §6.5 — 92% at 69–70 Mb vs 1% / 3% at 65–66 / 72–73 Mb |
| Which Trp53 haplotype survived (null vs wild-type) is not determined            | Needs client input | §7.2 — engineered "−" allele is coverage-invisible; needs its construct design |
| **Tumor3 has genome-wide copy-neutral LOH (haploidisation + endoreduplication)** |    High    | §6.6 — 95% LOH genome-wide, VAF bimodal at 0/1, copy number still ~2 |
| Study A tumors are aneuploid with distinct karyotypes                           |    High    | §6.2 |
| Study B tissues have stable diploid genomes; no aneuploidy at any age/genotype  |    High    | §6.2 — flat CN 1.94–2.08 |
| Oviduct phenotype is not driven by large-scale genomic instability              |    High    | §6.2 + validated positive control |
| iHPV transgene present specifically in L1L2H                                    |    High    | §6.5 |
| Study B background is not congenically pure C57BL/6J (~5–6M variants vs GRCm39) |    High    | §6.9 (129-type is the leading, unconfirmed, explanation) |
| Guide IDs are not GRCm39 coordinates (~76.5 kb offset)                          | High that they are not GRCm39; the mm10 attribution is inferred | §3.1 — no impact on results |

### 7.1 Implication for the model

**Tumor3 formed, and became the most aneuploid of the three tumors, while carrying none of the intended Brca1/Brca2/Pten edits.** The practical consequence is concrete: editing-escaped cells in the injected population can still form tumors, so a phenotype observed in a tumor should not be attributed to Brca/Pten loss without genotyping that specific tumor — as tumor3 demonstrates.

### 7.2 What drove the tumors: Trp53 loss, and in tumor3 a whole-genome catastrophe

**Every tumor lost Trp53 heterozygosity; none of the pre-tumor cells did** (§6.5). In tumor1 the homozygous tract is *focal* — 92% at 69–70 Mb where Trp53 sits, against 1% and 3% just 3 Mb to either side — which is selection for Trp53 loss written directly into the data. Tumor2 reaches the same endpoint through a large chr11 tract, tumor3 through genome-wide LOH. Three independent tumors, three different mechanisms, one converging target.

This supplies the driver the edits do not explain. Tumor3 carries **none** of the intended Brca1/Brca2/Pten lesions yet is the most aneuploid genome in the cohort; the resolution is that it took a different route to the same place — **whole-genome copy-neutral LOH** (§6.6), which renders the cell homozygous at every tumour suppressor simultaneously, Trp53 included. Tumors 1 and 2 combine their engineered Brca1+Pten knockouts with a Trp53 second hit; tumor3 dispensed with the engineering and lost Trp53 — and everything else — wholesale.

One limit, stated precisely: LOH establishes that each tumor retains a **single** Trp53 haplotype and is therefore homozygous, but *which* haplotype survived — the engineered null or the wild-type — cannot be read from these data, because the engineered lesion is too small to see in coverage (Trp53 read depth is flat in all six samples, including the parent). Retaining the wild-type copy would restore p53 function and confer no advantage, so loss of the wild-type allele is the overwhelmingly likely reading. If you send the construct design for that "−" allele, we can settle it directly on the existing data — no new sequencing needed.

### 7.3 Remaining analyses (no client input required unless noted)

1. **iHPV integration locus — the highest-value outstanding analysis for Study B.** Construct presence is confirmed, but where it integrated, and whether it disrupts a gene, is unresolved (§6.10, hypothesis 1). Needs the **full vector map** (PMC4662542) to capture junction-spanning reads. Note this can only ever explain L1L2H, not L1L2.
2. **Integrity of the lox-stop-lox cassette** (§6.10, hypothesis 2) — tests whether E6/E7 has been de-repressed by leaky recombination, and whether the Lats1/2 floxed exons are intact. Needs the **actual EGFP sequence** from the vector map; our marker used the wrong codon variant, so the assay has not yet run (§6.7).
3. **Lats1/2 expression in the affected tissue** — not a WGS analysis, but the one test that addresses the hypothesis most consistent with the data (a hypomorphic floxed allele affecting both strains, §6.10). Needs RNA or protein from oviduct.
4. **A wild-type littermate of the same background** — would separate "strain characteristic" from "consequence of the engineering" for the phenotype (§6.10). One WGS or even genotyping would do.
5. **Somatic point-mutation burden** — requires the additional filtering described in §6.11 before numbers are quoted.
6. **Study B de-novo candidate mining** — functional high-impact annotation (VEP/snpEff) + recurrence / genotype-differential filtering on the provided private call sets (§6.9). Recommended next analysis.
7. **Study B strain assignment** — optional; would settle whether the background is pure 129 or a 129×B6 mix.
8. **Which Trp53 haplotype the tumors retained** — the tumors are already shown to be homozygous at Trp53 (§6.5); confirming that the surviving copy is the engineered null rather than the wild-type needs the construct design of that "−" allele from you (§7.2). This is the only item requiring client input, and it refines an interpretation rather than changing a conclusion.

---

## 8. Deliverable Files

```
custom_research_report_20260716/
├── GeneEdit_Lats12_WGS_0716.md       ← this report
├── qc/                                MultiQC reports (Study A, Study B) — per-sample QC detail behind §6.1
├── cnv_ploidy/                        copy-number results behind §6.2
│   ├── cohort_cn_heatmap.png            all 12 samples, all chromosomes — one figure
│   ├── cohort_chrom_cn.tsv              chromosome-median copy number, per sample
│   ├── aneuploidy_calls.tsv             the gain/loss calls (CN ≥2.5 / ≤1.5)
│   └── *.cn_profile.png                 per-sample genome-wide CN profile (12 files)
├── edit_verification/                 guides, cut sites, cut-site genotypes, Brca2 allele quantification
│   ├── sgRNA_guides_reference.md        all 9 guides + Brca2 localisation / PAM check
│   ├── client_brca2_sgRNA_source_image.png   client image the Brca2 guides were transcribed from
│   ├── cut_sites.tsv                    all 9 cut sites (Pten/Brca1/Brca2) on GRCm39
│   ├── cutsite_indels.tsv               joint genotypes at every cut site, 6 samples
│   ├── brca2_allele_quant.tsv           per-sample wild-type / 31 bp / clip / indel counts
│   ├── brca2_clip_breakpoints.tsv       clip breakpoint stacking, per base per sample
│   ├── brca2_indels.tsv                 Brca2 gene-wide indel scan
│   └── spacer_gene.tsv                  spacer → target gene mapping
├── loh_trp53/                         loss-of-heterozygosity results (§6.5, §6.6)
│   ├── loh_analysis.txt                 genome-wide + per-chromosome LOH, VAF spectrum,
│   │                                    chr11 scan across Trp53 — with a how-to-read note
│   ├── chr11_scan_windows.bed           the chr11 windows scanned
│   └── trp53_coverage.tsv               Trp53 read depth (flat in all 6 — why LOH, not coverage,
│                                        is the informative assay)
├── identity_fingerprint/              tumor3 vs parent identity check
│   ├── fingerprint_summary.tsv          private / homozygous-private counts per sample
│   ├── fingerprint_breakdown.txt        per-chromosome + per-VAF breakdown & verdict
│   └── regions.bed                      the 38 × 500 kb sampled windows
├── ihpv_integration/                  construct-presence table, per-sample HPV16/EGFP/Luc reads (§6.5)
├── sv/sv_counts.tsv                   TIDDIT SV counts (§6.8) — raw per-sample, not somatic
├── somatic/
│   ├── somatic_counts.tsv               Mutect2 PASS counts — ⚠ artefact-inflated, NOT a burden (§6.10)
│   └── trp53_depth.tsv                  Trp53 locus depth — flat in all 6 samples; Trp53 status is
│                                        read from LOH instead (see loh_trp53/, §6.5)
└── candidates_denovo/                 MGP-filtered private variants (§6.9) — input for the recommended
                                        functional-annotation step, not a finished candidate list
```

This report and folder are complete and self-contained; they supersede the 2026-07-15 delivery, which is retained unchanged as the record of that round.

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics.*

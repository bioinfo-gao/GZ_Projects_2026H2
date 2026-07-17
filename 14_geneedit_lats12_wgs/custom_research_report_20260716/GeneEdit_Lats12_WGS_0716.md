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

1. **All intended CRISPR edits are verified, and the tumor lineages are resolved.** Brca1 + Pten knockout confirmed in B1TP. **Brca2 knockout confirmed in B2TP, and it is biallelic** — no wild-type Brca2 allele survives — alongside its Pten knockout. Tumors 1 and 2 derive from the **B1TP (Brca1+Pten)** lineage.
2. **Tumor 3 carries none of the intended edits and arose from an un-edited subclone of the parent.** It is wild-type at all three targeted genes (Brca1, Brca2, Pten), yet a SNP-fingerprint check shows it is genetically the **same animal** as the parental line — so it is an editing-escaped subclone, not a mis-tracked sample. **Practical consequence: editing-escaped cells can still form tumors in this model** (tumor3 is in fact the most aneuploid of the three), so a tumor phenotype should not be attributed to Brca/Pten loss without genotyping that specific tumor.
3. **Study A tumors are aneuploid, each with a distinct karyotype.** The three tumors carry different chromosome-scale gains/losses, consistent with independent clonal evolution. The parental and edited-but-pre-tumor cells remain diploid.
4. **Study B tissues are genomically stable — no aneuploidy.** All six Lats1/2-flox tissues (both genotypes, all ages) show a flat, diploid copy-number profile genome-wide (autosomal copy number 1.94–2.08). The reported oviduct phenotype is **not** driven by large-scale chromosomal instability or aneuploidy.
5. **The iHPV transgene is present specifically in L1L2H mice.** HPV16 E6/E7 and luciferase reads are detected in all three L1L2H samples and are entirely absent from all three L1L2 samples (perfect specificity), at levels consistent with a fixed germline transgene.
6. **Study B is not congenically pure C57BL/6J.** Each tissue carries ~5–6 million variants against GRCm39 — a genome-wide non-6J inbred background, most likely 129-derived from the ES-cell engineering. This is inherited background, not accumulated mutation, and it shapes how de-novo candidates must be mined (§6.6).

---

## 3. Sample Information

12 samples, paired-end 150 bp WGS (NovaSeq X Plus), aligned to GRCm39 (GENCODE vM35) with nf-core/sarek 3.8.1 (bwa-mem2). Mean autosomal depth in parentheses.

| #  | Label     | Study |  Type  | Group  | Genotype / role                                      |
| :- | :-------- | :---: | :----: | :----: | :--------------------------------------------------- |
| 1  | RO_origin |   A   |  Cell  | parent | Trp53⁺/⁻; Cas9; unedited — **matched normal** (~20×) |
| 2  | RO_B1TP   |   A   |  Cell  | edited | Brca1 + Pten KO (~24×)                               |
| 3  | RO_B2TP   |   A   |  Cell  | edited | Brca2 + Pten KO (~24×)                               |
| 4  | RO_tumor1 |   A   |  Cell  | tumor  | tumor from B1TP/B2TP injection (~26×)                |
| 5  | RO_tumor2 |   A   |  Cell  | tumor  | tumor (~19×)                                         |
| 6  | RO_tumor3 |   A   |  Cell  | tumor  | tumor (~30×)                                         |
| 7  | L1L2_3M   |   B   | Tissue |  L1L2  | Lats1/2 flox, 3 months (~21×)                        |
| 8  | L1L2H_3M  |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 3 months (~20×)                 |
| 9  | L1L2_12M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 12 months (~32×)                       |
| 10 | L1L2H_12M |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 12 months (~28×)                |
| 11 | L1L2_18M  |   B   | Tissue |  L1L2  | Lats1/2 flox, 18 months (~25×)                       |
| 12 | L1L2H_18M |   B   | Tissue | L1L2H  | Lats1/2 flox + iHPV, 18 months (~21×)                |

Total input: 12 samples, gzip FASTQ 462 GiB (measured on disk 2026-07-12).

### 3.1 CRISPR guides (CRISPRevolution sgRNA EZ Kit, 1.5 nmol each, Modified)

Nine guides supplied as images and converted RNA→DNA. Full table in `edit_verification/sgRNA_guides_reference.md`. The three Brca2 guides:

| Guide ID (as supplied) | Spacer (DNA, 5'→3')    | GRCm39 location              | Strand | PAM | Predicted cut    |
| :--------------------- | :--------------------- | :--------------------------: | :----: | :-: | :--------------: |
| Brca2+150529497        | `GATAAGCCTCAATTGGTTTG` | chr5:150,452,945–150,452,964 |   +    | AGG | chr5:150,452,961 |
| Brca2−150529492        | `AAAGCTCCTCAAACCAATTG` | chr5:150,452,954–150,452,973 |   −    | AGG | chr5:150,452,957 |
| Brca2−150529524        | `AGGTTCAGAATTGTATGGGG` | chr5:150,452,986–150,453,005 |   −    | GGG | chr5:150,452,989 |

Each maps to a single unique site inside Brca2, on the strand its name specifies, each with a perfect NGG PAM — confirming the guides are genuine and correctly transcribed. All three cut sites fall inside **Brca2 exon 3 (CDS)** of the Ensembl-canonical/CCDS transcript `ENSMUST00000044620.11` (*Brca2-201*, CCDS39411.1). Guides 1 and 2 overlap and cut the same point from opposite strands (4 bp apart); guide 3 cuts ~30 bp downstream — a multi-guide design concentrating three cuts into a **33 bp window**, which predicts either small indels or excision of the intervening fragment.

**Note on the guide IDs:** the coordinates embedded in the product names (e.g. `150529497`) **are not GRCm39 positions** — they sit ~76.5 kb from where the guides actually map on GRCm39, so using them directly against GRCm39 retrieves the wrong locus. The offset is most consistent with GRCm38/mm10 numbering, though this was not verified directly (no GRCm38 reference was used here). Every guide was located by sequence, so no result is affected.

---

## 4. Analysis Rationale and Decision Criteria

| Question                        | Approach & why it works without a dedicated wild-type control                                                                       | Threshold / criterion                                                      |
| :------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------- |
| Copy number / aneuploidy        | Ratio of genome-internal coverage → copy number; reference-free. GRCm39 = C57BL/6J ≈ the animals' own normal.                     | Chromosome-median CN ≥2.5 = gain, ≤1.5 = loss                            |
| Edit verification (A1)          | Read directly at each sgRNA cut site; compare edited/tumor vs the RO_origin parent. **Read/CIGAR level, not caller output** — see §6.2 | Wild-type reads absent ⇒ biallelic knockout; indel present in edited sample, absent (0/0) in RO_origin |
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
- **Edit verification:** multi-sample `bcftools mpileup`/`call` at cut-site windows (MAPQ≥20, BAQ≥15), **plus read-level allele quantification** (custom pysam analysis): reads with MAPQ≥20, duplicates excluded, classified by CIGAR into wild-type (spanning the cut-site core with ≥25 bp flanking margin, no indel/clip), the excision allele, other cut-site indels, or cut-site soft-clips (≥20 bp, breakpoint inside the core); clip breakpoints tabulated per base to separate recurrent from background clipping.
- **Sample identity fingerprint:** joint `bcftools mpileup/call` across all six Study A samples over 38 windows × 500 kb (19.0 Mb) spanning all 19 autosomes; biallelic SNPs retained. A "private allele" = sample carries ALT with ≥3 supporting reads while RO_origin has **0** ALT reads at ≥10× depth. Results broken down per chromosome and per allele fraction.
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

### 6.2 CRISPR edit verification (Study A)

**A methodological point that determines the result.** The joint `bcftools` call reported **no indel at any Brca2 cut site, in any sample**. Taken at face value this reads as "the Brca2 edit did not happen" — and it is wrong. At the B2TP cut site only **1 of 25** reads is a perfect match (`150M`), versus **24 of 34** in the unedited parent: the edited reads are present but carry a 31 bp deletion or are heavily soft-clipped (50–87 bp), which the caller did not resolve at this complex site. The Brca2 conclusions below therefore rest on read/CIGAR-level evidence.

Cut-site genotypes across the three target genes (RO_origin is homozygous reference throughout, confirming these are edits, not background):

| Gene / site             | RO_origin |     B1TP      |          B2TP          |         tumor1         |     tumor2      | tumor3 |
| :---------------------- | :-------: | :-----------: | :--------------------: | :--------------------: | :-------------: | :----: |
| Brca1 (chr11:101422906) |    0/0    | **indel**     |          0/0           | **indel (biallelic)**  | **indel (hom)** |  0/0   |
| Brca2 (chr5:150452957–989) |  wild-type |   wild-type   | **biallelic KO**  |       wild-type        |    wild-type    | wild-type |
| Pten (chr19:32777294)   |    0/0    | **indel**     | **indel (hom)**        |    **indel (hom)**     |      indel      |  0/0   |

**Brca2 allele composition at the cut site (read-level):**

| Sample        | Informative reads | **Wild-type** | 31 bp excision | Cut-site clips | Other indels             |
| :------------ | :---------------: | :-----------: | :------------: | :------------: | :----------------------- |
| RO_origin     |         7         |       6       |       0        |       1        | –                        |
| RO_B1TP       |        11         |      10       |       0        |       1        | –                        |
| **RO_B2TP**   |    **14**         |    **0**      |     **2**      |    **10**      | 3 bp del @150,452,989 ×2 |
| RO_tumor1     |        19         |      15       |       0        |       4        | –                        |
| RO_tumor2     |        12         |       8       |       0        |       4        | –                        |
| RO_tumor3     |        17         |       8       |       0        |       7        | 1 bp del @150,452,983 ×2 |

- **B1TP:** Brca1 KO + Pten KO — matches its Brca1+Pten design.
- **B2TP — Brca2 knockout confirmed, biallelic.** It is the **only sample with no wild-type Brca2 allele** (0 of 14 informative reads). Its Brca2 coverage is normal across the gene (23.2×) yet collapses to **0.57 of its own baseline** precisely at the cut window (13.3×) — the expected consequence of edited reads failing to align cleanly; every other sample sits at 0.87–1.23. Two distinct disrupted alleles are resolved:
  - **Allele 1 — a 31 bp deletion at chr5:150,452,958–150,452,988** (reads `76M31D43M31S`, `32S76M31D42M`), removing `TGGTTTGAGGAGCTTTCCTCAGAAGCCCCCC`. Its 5' end abuts the cut sites of the overlapping guide pair (150,452,957/961) and its 3' end abuts the third guide's cut (150,452,989): the deletion did not land merely "near" the guides — **its endpoints are the guides' cut sites**, i.e. exactly the fragment predicted to be excised when all three guides cut. It lies wholly within CDS exon 3, beginning at CDS nucleotide 91 (**codon ~31 of 3,329**), and 31 bp is not a multiple of 3 → **frameshift** → premature termination within the first 1% of the protein. This is a null allele, not a hypomorph.
  - **Allele 2 — disrupted at the third guide's cut site.** B2TP's 10 cut-site clips are not background: their breakpoints **stack on a single base** (150,452,990 ×4; 150,452,989 ×3 = 7 reads on that cut), whereas every other sample shows only scattered singletons. The clipped sequence maps only ambiguously elsewhere (MAPQ 0, repetitive), so the allele's precise structure — a complex indel or insertion — is not resolvable with short reads; that it is disrupted at the cut site is nonetheless clear.
- **B2TP = Brca2 + Pten knockout** is therefore established by direct observation of the Brca2 edit, not by inference from its Pten/Brca1 pattern.
- **tumor1, tumor2:** carry both Brca1 and Pten indels → **B1TP (Brca1+Pten) lineage**; both retain wild-type Brca2, as expected.
- **tumor3:** wild-type at **all three** targeted genes. Its only Brca2 indel call — a 1 bp deletion at 150,452,983 (2 reads) — falls inside a **7-C homopolymer** (`AAG`**`CCCCCCC`**`ATACAATTCTG`, 150,452,983–150,452,989), the classic context for a sequencing artifact, and is discounted. The lineage logic is decisive: **B2TP has no wild-type Brca2 allele**, so no descendant of B2TP could regain one; tumor3 is ~47% wild-type at this locus and therefore **cannot descend from B2TP**. Pten is equally decisive in the other direction — **both** injected lines carry a Pten edit (B2TP homozygously), so descent from either line would leave a detectable Pten lesion, and tumor3 has none (Pten sites 17–18× reads, unambiguous reference — a true absence of edit, not a coverage gap). Tumor3 descends from neither edited line.

### 6.3 Tumor3 identity — same animal as the parent

Tumor3 being wild-type at all three targeted genes leaves a question of **identity**: an un-edited escaper subclone of the parent, or a mis-tracked sample from a different mouse? B1TP/B2TP/tumor1/tumor2 — all known to descend from origin — provide the "same animal" baseline. 19,032 evaluable sites.

| Sample        | Private alleles (raw) | Excluding chr3/chr19 | **Homozygous private (excl.)** |
| :------------ | :-------------------: | :------------------: | :----------------------------: |
| RO_B1TP       |          245          |         223          |               0                |
| RO_B2TP       |          227          |         216          |               0                |
| RO_tumor1     |          264          |         231          |               0                |
| RO_tumor2     |          219          |         182          |               1                |
| **RO_tumor3** |      **2,851**        |         644          |            **8**               |

**The raw count is misleading, and is shown because it is instructive.** At face value tumor3's 2,851 private alleles (~11× every other sample) look like a different mouse. Two breakdowns show it is not:

- **By chromosome:** 1,154 of tumor3's 1,218 homozygous private alleles fall in a **single 500 kb window on chr3 (~47.9–48.4 Mb)** and 56 on chr19; **all 17 remaining autosomes contribute zero**. A genuinely different animal differs on *every* chromosome. This is a local event — and that chr3 window shows anomalous coverage in tumor3 alone (0.35× its own baseline, vs 0.69–1.00 in every other sample), i.e. a tumor-specific structural loss whose residual mismapped reads generate spurious homozygous calls.
- **By allele fraction:** excluding those two regions, tumor3 retains just **8 homozygous private alleles**. Its remaining excess sits at **VAF 0.3–0.7** — the signature of **clonal somatic mutation**, expected in the most rearranged of the three tumors, not of inherited variation, which would pile up at VAF ~1.0.

**Conclusion:** tumor3 shares RO_origin's germline genome. It is the same animal, and a swap with a different mouse is excluded. With §6.2, **tumor3 arose from an un-edited (editing-escaped) subclone of the Trp53⁺/⁻;Cas9 parent** — and its heavy aneuploidy therefore developed independently of any Brca/Pten editing.

**One caveat, stated plainly:** this method distinguishes *different genomes*, not *different individuals of an inbred line*. Littermates of an inbred colony are near-genetically identical, so a mix-up between two such animals would be invisible to any fingerprint. What is established is that tumor3's genome is the parent's genome — precisely what the un-edited-subclone explanation requires, and what a swap with an unrelated or differently-engineered animal would have broken.

### 6.4 iHPV transgene detection (Study B)

Construct-marker reads (unmapped reads aligned to HPV16 + EGFP + luciferase):

| Sample    | Construct reads | HPV16 | Luciferase |
| :-------- | :-------------: | :---: | :--------: |
| L1L2_3M   |        0        |   0   |     0      |
| L1L2H_3M  |       45        |  16   |     71     |
| L1L2_12M  |        0        |   0   |     0      |
| L1L2H_12M |       83        |  41   |    123     |
| L1L2_18M  |        0        |   0   |     0      |
| L1L2H_18M |       58        |  37   |     75     |

**The iHPV transgene is present in all three L1L2H samples and absent in all three L1L2 samples** — a perfectly specific result. Read counts do not increase with age, consistent with a fixed germline transgene rather than an age-accumulating event.

### 6.5 Structural variants

TIDDIT SV counts are available for all samples (`sv/sv_counts.tsv`). For Study B these are dominated by C57BL/6-versus-reference background SVs and are being refined with the MGP background filter before interpretation; Study A tumor SVs will be reported as tumor-minus-normal somatic events.

### 6.6 De-novo candidate variants (Study B)

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

Genuine causal de-novo events number in the dozens-to-hundreds, not hundreds of thousands. **Identifying them from this background requires functional restriction to high-impact coding consequences (frameshift / stop-gain / splice, via VEP or snpEff) combined with recurrence and L1L2-vs-L1L2H differential filtering.** That functional-annotation step is the recommended immediate next analysis; the per-sample private call sets are provided for it (`candidates_denovo/*.private.vcf.gz`). Critically, this does not affect the primary Study B conclusion: the copy-number analysis (§6.1) already establishes these genomes are structurally stable and diploid.

### 6.7 Somatic point mutations (Study A)

Mutect2 tumor-vs-origin calling was run for all five pairs. **The current PASS counts are inflated by artifacts (no contamination table / panel-of-normals) and are not yet a reliable mutation burden;** they are therefore not reported here as findings. This analysis requires additional filtering (VAF/depth thresholds, orientation-bias, PoN) before mutation burdens are quoted. The Study A tumor genome story is presently carried by the copy-number (§6.1), edit-verification (§6.2) and identity (§6.3) results, which are robust.

---

## 7. Conclusions

| Conclusion                                                                     | Confidence | Evidence |
| :----------------------------------------------------------------------------- | :--------: | :------- |
| The three Brca2 guides are genuine and target Brca2 exon 3                     |    High    | §3.1 — unique GRCm39 hit each, correct strand, perfect NGG PAM |
| **B2TP carries the intended Brca2 knockout, and it is biallelic**              |    High    | §6.2 — 0/14 wild-type reads; coverage 0.57× of own baseline at the cut window |
| Brca2 allele 1 = 31 bp dual-cut excision → frameshift at codon ~31/3,329 → null |    High    | §6.2 — deletion endpoints coincide with guide cut sites |
| Brca2 allele 2 = disrupted at the third guide's cut; exact structure unresolved | Moderate–High | §6.2 — 7-read breakpoint stack on 150,452,989/990 |
| **B2TP = Brca2 + Pten knockout**, by direct observation rather than inference   |    High    | §6.2 |
| B1TP = Brca1 + Pten KO confirmed                                                |    High    | §6.2 |
| Tumors 1 & 2 arose from the B1TP lineage; neither carries a Brca2 edit          |    High    | §6.2 |
| **Tumor3 arose from an un-edited subclone of the parent, not B1TP/B2TP**        |    High    | §6.2 (wild-type at all three genes) + §6.3 (shares origin's germline) |
| **Tumor3 is the same animal as RO_origin — sample swap excluded**               | High (see §6.3 caveat re inbred littermates) | §6.3 — 8 homozygous private alleles vs thousands expected for a different mouse |
| Tumor3's chr3 ~47.9–48.4 Mb anomaly is a tumor-specific structural event, not an identity signal | Moderate–High | §6.3 — confined to 1 window; coverage 0.35× of own baseline |
| Study A tumors are aneuploid with distinct karyotypes                           |    High    | §6.1 |
| Study B tissues have stable diploid genomes; no aneuploidy at any age/genotype  |    High    | §6.1 — flat CN 1.94–2.08 |
| Oviduct phenotype is not driven by large-scale genomic instability              |    High    | §6.1 + validated positive control |
| iHPV transgene present specifically in L1L2H                                    |    High    | §6.4 |
| Study B background is not congenically pure C57BL/6J (~5–6M variants vs GRCm39) |    High    | §6.6 (129-type is the leading, unconfirmed, explanation) |
| Guide IDs are not GRCm39 coordinates (~76.5 kb offset)                          | High that they are not GRCm39; the mm10 attribution is inferred | §3.1 — no impact on results |

### 7.1 Implication for the model

**Tumor3 formed, and became the most aneuploid of the three tumors, while carrying none of the intended Brca1/Brca2/Pten edits.** The practical consequence is concrete: editing-escaped cells in the injected population can still form tumors, so a phenotype observed in a tumor should not be attributed to Brca/Pten loss without genotyping that specific tumor — as tumor3 demonstrates.

### 7.2 What drove tumor3's aneuploidy is not established

Trp53 loss is the obvious candidate, given the Trp53⁺/⁻ parental genotype. **We could not test it.** A copy-number check across Trp53 (chr11:69,471,185–69,482,699) returns a flat, unchanged profile in **all six samples, including the parent**. Since RO_origin is Trp53⁺/⁻ by design yet shows no coverage deficit, its engineered null allele must be a small, coverage-invisible lesion rather than a deletion — so this assay cannot read Trp53 status in any sample, and no conclusion is drawn from it. To determine whether tumor3 lost its remaining Trp53 allele we would need the design of that "−" allele (what the engineered lesion actually is); with it, the existing data are sufficient to genotype it directly.

### 7.3 Remaining analyses (no client input required unless noted)

1. **iHPV integration locus** — construct presence is confirmed, but the base-pair integration junction / disrupted gene is not yet resolved; this requires the full PMC4662542 vector map to capture junction-spanning reads.
2. **Somatic point-mutation burden** — requires the additional filtering described in §6.7 before numbers are quoted.
3. **Study B de-novo candidate mining** — functional high-impact annotation (VEP/snpEff) + recurrence / genotype-differential filtering on the provided private call sets (§6.6). Recommended next analysis.
4. **Study B strain assignment** — optional; would settle whether the background is pure 129 or a 129×B6 mix.
5. **Trp53 genotyping in tumor3** — needs the Trp53 "−" allele design from you (§7.2).

---

## 8. Deliverable Files

```
custom_research_report_20260716/
├── GeneEdit_Lats12_WGS_0716.md       ← this report
├── qc/                                MultiQC reports (Study A, Study B)
├── cnv_ploidy/                        cohort CN table, aneuploidy calls, heatmap, per-sample CN profiles
├── edit_verification/                 guides, cut sites, cut-site genotypes, Brca2 allele quantification
│   ├── sgRNA_guides_reference.md        all 9 guides + Brca2 localisation / PAM check
│   ├── client_brca2_sgRNA_source_image.png   client image the Brca2 guides were transcribed from
│   ├── cut_sites.tsv                    all 9 cut sites (Pten/Brca1/Brca2) on GRCm39
│   ├── cutsite_indels.tsv               joint genotypes at every cut site, 6 samples
│   ├── brca2_allele_quant.tsv           per-sample wild-type / 31 bp / clip / indel counts
│   ├── brca2_clip_breakpoints.tsv       clip breakpoint stacking, per base per sample
│   ├── brca2_indels.tsv                 Brca2 gene-wide indel scan
│   └── spacer_gene.tsv                  spacer → target gene mapping
├── identity_fingerprint/              tumor3 vs parent identity check
│   ├── fingerprint_summary.tsv          private / homozygous-private counts per sample
│   ├── fingerprint_breakdown.txt        per-chromosome + per-VAF breakdown & verdict
│   └── regions.bed                      the 38 × 500 kb sampled windows
├── ihpv_integration/                  construct-presence table (per-sample HPV16/EGFP/Luc reads)
├── somatic/                           Mutect2 PASS counts (preliminary), Trp53 locus depth
├── sv/                                TIDDIT SV counts
└── candidates_denovo/                 MGP-filtered de-novo candidates (for the §6.6 next step)
```

This report and folder are complete and self-contained; they supersede the 2026-07-15 delivery, which is retained unchanged as the record of that round.

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics.*

# Shotgun Metagenomics of Mouse Gut Microbiome under Ad Libitum vs Intermittent Fasting (High-Fat Diet)

**Project:** QTE_26_06_25_001_Daniel_Mendes
**Report Date:** 2026-07-18
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Mus musculus* (host genome GRCm39)
**Tissue / Material:** Stool (fecal pellets)
**Sequencing:** Illumina NovaSeq X Plus, paired-end 150 bp, shotgun metagenomics

> **Scope of this report.** This document covers the **assembly-free (reference-based) analysis**: taxonomic composition, community diversity, and **functional pathway profiling (HUMAnN)** — the full standard shotgun-metagenomics deliverable. A genome-resolved **metagenome-assembled genome (MAG)** analysis (Phase 2) is running as a separate, genome-level extension and will be delivered subsequently.

---

## 1. Objectives

1. Characterize the fecal microbiome composition of high-fat-diet (HFD) mice under two feeding regimens — **ad libitum (AL)** vs **intermittent fasting (IF)** — at genus and species resolution.
2. Compare **alpha diversity** (within-sample richness/evenness) and **beta diversity** (between-sample community structure) between the two diet arms.
3. Identify bacterial taxa whose relative abundance differs between AL and IF.
4. Cross-validate the taxonomic profile using two independent classifiers (Kraken2/Bracken and MetaPhlAn).

## 2. Key Findings

- **Deep, high-quality data.** 227.5 million read pairs total (~68 Gbp; 17–27 M pairs/sample). Host (mouse) contamination was low (2.2–18.2%, mean ~6.8%), as expected for stool, leaving abundant microbial signal.
- **The two diet arms are not separable by overall community structure.** Bray-Curtis PERMANOVA was not significant (R² = 0.149, *p* = 0.25); between-group community distance (0.295) essentially equals within-group distance (0.298). With n = 5 per arm and high inter-individual variability, HFD is the dominant driver and the IF vs AL effect on global structure is not statistically resolvable in this cohort.
- ***Akkermansia muciniphila* dominates and trends higher under IF.** It is the single most abundant species (mean 48.7% in AL vs 64.4% in IF). A higher *Akkermansia* abundance under fasting is consistent with the published intermittent-fasting literature, but here the difference is **not statistically significant** (Wilcoxon *p* = 0.42).
- **No taxon reaches statistical significance after multiple-testing correction.** Among 60 core abundant species, none passed FDR < 0.05. Suggestive (non-significant) trends: *Lactococcus* / *Lactobacillus johnsonii* higher in AL; *Akkermansia* and *Paramuribaculum* higher in IF.
- **Alpha diversity trends slightly higher in AL** across all three indices (Observed, Shannon, Simpson) but not significantly (*p* = 0.31–0.42). IF mice were more homogeneous among themselves (within-group distance 0.256 vs 0.340 for AL).
- **Functional (metabolic pathway) profiles agree with the taxonomic picture: no *statistically resolvable* diet effect (underpowered, not a demonstrated null).** HUMAnN pathway profiles do not separate AL from IF (PERMANOVA R² = 0.155 — diet still accounts for ~16% of functional variance — *p* = 0.22), and no pathway is significant after FDR correction. The consistent (non-significant) trend is a mild shift toward higher biosynthetic capacity under IF (amino-acid, nucleotide, and peptidoglycan biosynthesis) versus slightly higher glycolysis under AL.

## 3. Sample Information

Ten fecal samples, five per diet arm. Replicate identifiers (e.g. `4_02_25`, `6_05_12`) recur across both arms; these are treated here as **independent biological replicates** (two-group design). If these identifiers denote litter/cage/batch matching, a paired analysis can be run on request (see §5).

| Sample | Diet arm | Raw read pairs | Host (mouse) reads removed |
| :--- | :---: | :---: | :---: |
| HFD_AL_4_02_25 | AL | 24.0 M | 6.35% |
| HFD_AL_4_03_11 | AL | 22.9 M | 5.23% |
| HFD_AL_6_05_12 | AL | 24.4 M | 18.19% |
| HFD_AL_6_05_22 | AL | 26.5 M | 10.06% |
| HFD_AL_7_06_12 | AL | 19.1 M | 6.07% |
| HFD_IF_4_02_25 | IF | 21.2 M | 7.76% |
| HFD_IF_4_03_11 | IF | 17.2 M | 2.20% |
| HFD_IF_6_05_12 | IF | 24.0 M | 2.31% |
| HFD_IF_6_05_22 | IF | 24.6 M | 6.04% |
| HFD_IF_7_06_12 | IF | 23.7 M | 3.40% |

Total: **227.5 M read pairs (~68 Gbp)**.

## 4. Analysis Rationale and Decision Criteria

- **Assembly-free, reference-based route chosen as the primary analysis.** The scientific question is a case-control community comparison. Mouse gut is well represented in reference databases, so read-level classification answers composition, diversity, and (forthcoming) function directly and robustly, without the depth demands and interpretive overhead of de-novo assembly.
- **Two independent classifiers for cross-validation.** *Kraken2 + Bracken* (k-mer based, whole-community, casts a wide net) and *MetaPhlAn 4* (clade-specific marker genes, conservative). Concordant calls are trustworthy; discordances are flagged rather than hidden.
- **Host removal against GRCm39** before classification, so mouse reads do not inflate microbial abundance. Residual human reads (handling/kit contamination) were additionally removed from the abundance analysis.
- **Differential-abundance testing restricted to a "core" set** — species with mean relative abundance ≥ 0.1% detected in ≥ 5/10 samples. This is a deliberate decision: k-mer classifiers assign a long tail of spurious near-zero taxa (here 5,479 raw Bracken species vs 471 confidently detected by MetaPhlAn), and testing that tail produces large but meaningless fold-changes on environmental/contaminant organisms. Core-filtering keeps the biology and removes the noise.
- **Significance threshold:** Wilcoxon rank-sum test, Benjamini-Hochberg FDR < 0.05. Beta diversity assessed by PERMANOVA (999 permutations) on Bray-Curtis distances.

## 5. Methods

| Step | Tool / parameters |
| :--- | :--- |
| Read QC & trimming | fastp (adapter/quality trimming) |
| Host removal | Bowtie2 vs *M. musculus* GRCm39 primary assembly; unmapped (microbial) reads retained |
| Taxonomic profiling | nf-core/taxprofiler 2.0.1; **Kraken2** + **Bracken** (`-r 150`) on Standard-8GB DB; **MetaPhlAn 4** (CHOCOPhlAn SGB vJan25) |
| Profile merging | taxpasta; MetaPhlAn/Bracken combined abundance tables |
| Diversity | vegan (Shannon, Simpson, observed richness; Bray-Curtis; PCoA; adonis2 PERMANOVA) |
| Differential abundance | Wilcoxon rank-sum on relative abundance, Benjamini-Hochberg FDR; core filter (≥0.1% mean, ≥5/10 prevalence) |
| Functional profiling | HUMAnN 3.9 (ChocoPhlAn + UniRef90, vJun23 taxonomic prescreen); MetaCyc pathway & UniRef90 gene-family abundance; community-level pathways compared AL vs IF (Wilcoxon + BH, core filter; Bray-Curtis PCoA + PERMANOVA) |

Analysis was orchestrated with resource caps within server policy (≤ 28 cores). All figures are provided as both 300-dpi PNG and vector PDF.

## 6. Results

### 6.1 Data quality and host content
All ten libraries were deep (17–27 M read pairs) and of high quality. Host (mouse) contamination was low and variable (2.2–18.2%), consistent with fecal material; the majority of reads in every sample were microbial and available for classification.

### 6.2 Community composition (Figure 1)
The mouse fecal community is dominated by *Akkermansia*, unclassified/uncultured genome bins (MetaPhlAn SGB placeholders labeled `GGB…`), *Lachnospiraceae*, and *Oscillibacter*, with *Dubosiella*, *Muribaculum*, *Faecalibaculum*, and *Duncaniella* — a typical HFD mouse gut profile. Composition is broadly similar between arms, with substantial mouse-to-mouse variation within each arm (e.g. *Akkermansia* ranges widely across individuals). *Akkermansia muciniphila* is the single most abundant species (mean 48.7% AL, 64.4% IF by Bracken).

### 6.3 Alpha diversity (Figure 2)
Ad-libitum samples trended toward higher within-sample diversity than intermittent-fasting samples across all three indices (Observed richness, Shannon, Simpson), but none of the differences were statistically significant (Wilcoxon *p* = 0.31, 0.31, 0.42 respectively). Note that observed-richness values are inflated by the k-mer classifier's low-abundance tail and should be read as a relative comparison, not an absolute species count; the abundance-weighted Shannon and Simpson indices are more robust and show the same non-significant trend.

### 6.4 Beta diversity (Figure 3) — *ordination interpretation*
Bray-Curtis PCoA places the first two axes at **66.9% (PCo1)** and **16.8% (PCo2)** of variance (83.7% combined). **The AL and IF samples do not separate:** their 80% confidence ellipses overlap almost entirely, and quantitatively the mean between-arm community distance (0.295) is indistinguishable from the mean within-arm distance (0.298) — a between/within ratio of 0.99. PERMANOVA finds no statistically reliable partitioning by diet (**R² = 0.149, *p* = 0.25**, 999 permutations); diet explains ~15% of community variation — a modest but non-zero effect that is not statistically resolvable at this sample size (an underpowered test, not evidence of no effect; see §7.1). One secondary observation: the IF mice are more similar to one another (within-group distance 0.256) than the AL mice are (0.340), i.e. fasting may tighten the community toward a more consistent state, though this is descriptive only.

### 6.5 Differential abundance (Figure 4)
Across the 60 core abundant species, **no species reached FDR < 0.05** (smallest adjusted *p* = 0.60). The strongest *trends* (all non-significant) were: *Akkermansia muciniphila* and *Paramuribaculum intestinale* higher under IF; *Lactobacillus johnsonii*, *Lactococcus* spp., and *Clostridium* sp. higher under AL. The *Akkermansia* trend is directionally consistent with the fasting/metabolic-health literature and is the most biologically noteworthy signal, but the present cohort (n = 5/arm, high individual variability) is not powered to establish it statistically.

### 6.6 Cross-tool concordance (Figure 5)
Bracken and MetaPhlAn agreed moderately on species-level mean abundances (Spearman ρ = 0.64). Bracken systematically reported higher abundances than MetaPhlAn, and the two tools disagreed on the identity of the dominant taxon: Bracken assigns ~50% to *Akkermansia muciniphila*, whereas MetaPhlAn distributes much of the dominant signal to unnamed species-level genome bins (SGBs). This is a known consequence of the two databases' differing taxonomies and should be kept in mind when quoting absolute percentages; both tools agree the community is *Akkermansia*-rich.

### 6.7 Functional pathway profiling (Figures 6–8)
HUMAnN reconstructed **1,243 MetaCyc pathway features** across the cohort (community-level and species-stratified); analysis focused on the 128 well-represented community pathways. The functional picture closely mirrors the taxonomic one:
- **Beta diversity (Figure 6):** pathway-level Bray-Curtis PCoA shows AL and IF overlapping heavily, with no statistically reliable separation (**PERMANOVA R² = 0.155, *p* = 0.22**; PCo1 = 59.2%). Note that R² = 0.155 means diet still accounts for ~16% of the between-sample variation in functional potential — a modest but non-negligible effect size that this 5-vs-5 design is underpowered to resolve; the non-significance reflects limited power, not a demonstrated absence of effect (see §7.1, Recommendations).
- **Differential pathways (Figure 7):** **no pathway is significant after FDR correction** (smallest adjusted *p* = 0.76). The directional trends are biologically coherent: biosynthetic pathways — UDP-N-acetylmuramoyl-pentapeptide / peptidoglycan biosynthesis, branched-chain and aromatic amino-acid biosynthesis, UMP and purine nucleotide biosynthesis — trend modestly higher under **IF**, while glycolysis IV and one arginine-biosynthesis route trend higher under **AL**. This is consistent with a mild fasting-associated shift toward biosynthetic/anabolic capacity, but it is a trend only, not a significant finding.
- **Pathway landscape (Figure 8):** the top-25 most abundant pathways (dominated by nucleotide, amino-acid, and cofactor biosynthesis and central carbon metabolism) are present and comparably abundant in both arms, confirming a shared core metabolic repertoire.

## 7. Conclusions

| Question | Result |
| :--- | :--- |
| Overall community structure differs AL vs IF? | **No statistically reliable difference** — PERMANOVA n.s. (R² = 0.149, i.e. diet accounts for ~15% of variance, *p* = 0.25); between ≈ within-group distance. This is an *underpowered* test at n = 5/arm, not a demonstrated absence of effect |
| Alpha diversity differs? | **No statistically reliable difference**; AL trends slightly higher and in the *same direction across all three indices* (*p* = 0.31–0.42) — a consistent but underpowered signal |
| Any taxon significantly differential? | **None reaches FDR < 0.05** (n = 5/arm is underpowered for per-taxon testing after multiple-testing correction); directional trends only |
| Dominant taxon | *Akkermansia muciniphila* (~49% AL, ~64% IF), trending higher under IF (n.s.) |
| Notable non-significant trends | ↑IF: *Akkermansia*, *Paramuribaculum*; ↑AL: *Lactobacillus johnsonii*, *Lactococcus* |
| Cross-tool agreement | Moderate (Spearman ρ=0.64); both agree community is *Akkermansia*-rich |
| Functional pathways differ AL vs IF? | **No statistically reliable difference** — PERMANOVA n.s. (R² = 0.155, i.e. diet accounts for ~16% of functional variance, *p* = 0.22); no pathway at FDR < 0.05. Underpowered, not a demonstrated null |
| Functional trend | Mild ↑ biosynthesis (amino acid/nucleotide/peptidoglycan) in IF; ↑ glycolysis in AL (n.s.) |

**Overall interpretation.** In this HFD cohort, intermittent fasting did not produce a *statistically resolvable* difference in the fecal microbiome relative to ad-libitum feeding — at the level of community structure, species composition, or metabolic-pathway potential. Importantly, this should be read as an **underpowered result, not as evidence that diet has no effect**: across every layer the estimated effect size is modest but non-zero and directionally coherent (community R² = 0.149, functional R² = 0.155 — i.e. diet accounts for ~15–16% of variance; the *Akkermansia*/biosynthesis trends all point the same way), and a 5-vs-5 design with high inter-individual variability is simply not powered to resolve an effect of this size. The dominant, biologically plausible signal is a fasting-associated enrichment of *Akkermansia muciniphila* with a coherent mild shift toward biosynthetic pathway capacity. We regard these as **real candidate effects worth confirming**, and outline below how a follow-up could be powered to test them. The genome-resolved MAG analysis (Phase 2, in progress) will separately test whether diet effects are more apparent at the strain/genome level than at the species-composition or pathway level.

### 7.1 Statistical power and recommendations for a follow-up study

The recurring theme across all analyses is **low statistical power**, not absence of biological effect. An effect explaining ~15–16% of community/functional variance (R² ≈ 0.15) is, when significant, a publishable microbiome effect; here it falls short of significance only because n = 5/arm against high individual variability leaves little power. A confirmatory study could be designed to resolve these trends by acting on the three levers that drive PERMANOVA / rank-test power:

1. **Increase replication.** This is the single most reliable lever. Roughly doubling to **≈10–12 animals per arm** would substantially raise power for an effect of the observed size; a formal power estimate can be produced from the current variance if a target effect size is agreed.
2. **Strengthen and standardize the intervention.** A longer or more strictly enforced fasting window widens the AL-vs-IF physiological contrast, enlarging the between-arm effect (R²) rather than only fighting variance. Standardizing the regimen across animals also shrinks within-arm spread.
3. **Reduce within-arm variance.** Co-house or use cage-/litter-matched littermates to control the strong **cage/coprophagy effect** (in mice the cage is often the single largest source of gut-microbiome variance), and collect stool at a **fixed point in the fasting cycle** so each animal's snapshot is comparable — this matters *more for the functional than the taxonomic layer*, because gene/pathway activity tracks feeding state. The IF arm is already the more homogeneous of the two (within-arm Bray–Curtis 0.256 vs 0.340 for AL), so tightening the **AL arm** in particular should yield the most gain.

**Leverage the existing data at no extra cost:** if the recurring replicate identifiers (e.g. `4_02_25`, `6_05_12`) correspond to cage- or litter-matched pairs across the two arms, a **paired / blocked test** (PERMANOVA with the pair as a stratum, or a paired Wilcoxon) would already recover power from the samples in hand — we can re-run this immediately if the pairing is confirmed. Note also that even with a significant global test, individual pathway/taxon calls face multiple-testing correction and need their own power, so per-feature claims should still be interpreted cautiously.

## 8. Deliverable Files

```
custom_research_report_20260718/
├── Daniel_Mendes_gut_metagenomics_0718.md   ← this report
├── qc/
│   └── multiqc_report.html                  ← per-sample QC (fastp, FastQC, host removal, classifier stats)
├── taxonomy/
│   ├── fig1_composition_genus.(png|pdf)      ← genus stacked bars, by arm
│   ├── fig1b_composition_species.(png|pdf)   ← species stacked bars, by arm
│   ├── fig5_crosstool_concordance.(png|pdf)  ← Bracken vs MetaPhlAn
│   ├── composition_genus_metaphlan.tsv
│   ├── composition_species_metaphlan.tsv
│   ├── bracken_species_abundance.txt         ← full Bracken species table
│   └── metaphlan_abundance.txt               ← full MetaPhlAn profile
├── diversity/
│   ├── fig2_alpha_diversity.(png|pdf)
│   ├── fig3_beta_pcoa.(png|pdf)
│   ├── fig4_differential_abundance.(png|pdf)
│   ├── braycurtis_distance.tsv
│   ├── differential_abundance_CORE_bracken.tsv   ← core species tested
│   └── differential_abundance_ALL_bracken.tsv    ← full table (for reference)
└── function/                                 ← HUMAnN functional profiling
    ├── fig6_functional_pcoa.(png|pdf)        ← pathway Bray-Curtis PCoA + PERMANOVA
    ├── fig7_pathway_differential.(png|pdf)   ← diet-associated pathway trends
    ├── fig8_top_pathways_heatmap.(png|pdf)   ← top-25 pathway landscape
    ├── differential_pathways_humann.tsv      ← pathway Wilcoxon/BH results
    ├── pathway_abundance_relab.tsv           ← MetaCyc pathway relative abundance (all samples)
    └── genefamilies_cpm.tsv                  ← UniRef90 gene-family abundance (CPM)
```

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics — 2026-07-18*

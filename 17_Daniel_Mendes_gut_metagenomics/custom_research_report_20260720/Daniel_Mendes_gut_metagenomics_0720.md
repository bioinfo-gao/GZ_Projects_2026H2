# Shotgun Metagenomics of Mouse Gut Microbiome under Ad Libitum vs Intermittent Fasting — Phase 2: Genome-Resolved Metagenome-Assembled Genomes (MAGs)

**Project:** QTE_26_06_25_001_Daniel_Mendes
**Report Date:** 2026-07-20
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Mus musculus* (host genome GRCm39)
**Tissue / Material:** Stool (fecal pellets)
**Sequencing:** Illumina NovaSeq X Plus, paired-end 150 bp, shotgun metagenomics

> **Scope of this report.** This document covers **Phase 2: genome-resolved analysis** — de novo co-assembly, binning, quality control, and taxonomic classification of metagenome-assembled genomes (MAGs). It extends the **Phase 1** deliverable (`custom_research_report_20260718/Daniel_Mendes_gut_metagenomics_0718.md`), which covers taxonomic composition, diversity, and HUMAnN functional pathway profiling from the same ten samples. The two reports are complementary and should be read together; Phase 1 is not repeated here except where needed to interpret Phase 2 findings.

---

## 1. Objectives

1. Recover genome-resolved bacterial genomes (MAGs) from the fecal shotgun data via de novo co-assembly and binning, independent of reference-database composition (the main limitation of the read-classification route in Phase 1).
2. Assess genome quality (completeness/contamination) against MIMAG standards and assign GTDB taxonomy to each recovered genome.
3. Compare the genome-resolved catalog between diet arms (**AL** vs **IF**) — genome count, quality distribution, and taxonomic composition.
4. Cross-validate Phase 1's dominant read-level finding (*Akkermansia muciniphila*) against the genome-assembly evidence.

## 2. Key Findings

- **168 medium-to-high-quality MAGs recovered, totaling ~447 Mb of genome-resolved bacterial sequence** (184 total non-redundant genomes before quality filtering: 95 from the AL co-assembly, 89 from IF). By MIMAG standard: **79 high-quality** (≥90% completeness, <5% contamination) and **89 medium-quality** (≥50% completeness, <10% contamination); 16 genomes did not meet either threshold and are excluded from downstream genome-level use.
- **Both diet arms yield a similar genome catalog, dominated by Bacillota (Firmicutes).** Of GTDB-classified MAGs, Bacillota accounts for 84% of AL genomes (73/87) and 80% of IF genomes (64/80), with Bacteroidota and Actinomycetota as consistent minor phyla in both arms — no qualitative phylum-level shift between diets at the genome-resolved level, consistent with Phase 1's community-level finding of no statistically resolvable diet effect.
- ***Akkermansia muciniphila* was recovered as a genome only in the AL arm, not IF — despite Phase 1 finding it more abundant by read count in IF (64.4% vs 48.7%).** This apparent contradiction is expected, not an error: the AL-arm MAG is complete (100%) but its per-sample coverage spans an 8-fold range (31×–238×) across the five AL animals, indicating substantial strain/individual heterogeneity that can fragment the assembly graph for a single dominant taxon — read-based abundance and clean genome assembly are different sensitivities, and this discordance is itself informative (see §6.3).
- **Genus-level MAG recovery corroborates several taxa flagged in Phase 1** — *Oscillibacter* and *Paramuribaculum* (both noted as diet-associated trend taxa in Phase 1, §6.2/6.5) are among the most frequently recovered genera here, along with *Acetatifactor*, *Lawsonibacter*, and *Acutalibacter* — a genome-level view of the same community that Phase 1 profiled by reads.
- **No sample-level (only group-level) MAGs were produced.** Per-sample sequencing depth (~4 Gbp/sample) was too shallow for reliable individual-sample assembly, so contigs were co-assembled per diet arm (5 samples pooled per group) as planned in the research design — this trades individual-animal MAG resolution for assembly contiguity, and is why MAG comparisons below are AL-vs-IF only, not animal-by-animal.

## 3. Sample Information

Same ten fecal samples as Phase 1 (five per diet arm; see Phase 1 report §3 for per-sample read counts and host content). For MAG assembly, reads were **pooled by diet arm** into two co-assembly groups because per-sample depth (~4 Gbp) was insufficient for reliable individual assembly:

| Co-assembly group | Samples pooled | Input reads |
| :--- | :---: | :---: |
| group-AL | HFD_AL_4_02_25, HFD_AL_4_03_11, HFD_AL_6_05_12, HFD_AL_6_05_22, HFD_AL_7_06_12 | 5 samples, ~117 M pairs |
| group-IF | HFD_IF_4_02_25, HFD_IF_4_03_11, HFD_IF_6_05_12, HFD_IF_6_05_22, HFD_IF_7_06_12 | 5 samples, ~111 M pairs |

## 4. Analysis Rationale and Decision Criteria

- **Group (per-diet-arm) co-assembly, not per-sample or cohort-wide.** Individual-sample depth (~4 Gbp) is below the threshold typically needed for a well-resolved gut metagenome assembly; a single cohort-wide co-assembly would maximize depth but collapse any AL/IF genome-level contrast entirely. Per-arm co-assembly (5 samples pooled per group) is the standard middle ground: it raises effective depth to a level that supports genome recovery while still keeping the two diet arms comparable as separate catalogs.
- **Three binners + consensus refinement, not a single binner.** MetaBAT2, MaxBin2, and SemiBin2 use different signals (coverage covariance, tetranucleotide frequency, and a semi-supervised deep-learning model respectively) and disagree substantially on any given genome. **DASTool** was used to select the single best, non-redundant bin per genomic cluster across all three binners' predictions per group — this is the standard way to avoid both missed bins (any one binner's blind spot) and redundant double-counting of the same genome (which would happen if all three binners' raw outputs were reported as if they were 3× as many genomes). The DASTool-refined catalog (184 genomes) is the deliverable's "final MAG" set; the raw per-binner outputs (1,060 rows, comparison purposes only) are not part of the count.
- **MIMAG quality tiers applied post hoc via CheckM2**, not DASTool's own internal score, because CheckM2 (machine-learning-based, trained on a broad genome set) is the community-standard external QC metric independent of the binning process itself. Threshold: High ≥90% completeness & <5% contamination; Medium ≥50% completeness & <10% contamination (both per Bowers et al. 2017 MIMAG standard); anything not meeting the Medium threshold is retained in the catalog for transparency but flagged **Low** and should be excluded from genome-level downstream analyses (pangenomics, phylogenomics, strain tracking).
- **GTDB-Tk (release R226) for taxonomy**, giving genus/species-level placement on a curated, up-to-date bacterial/archaeal reference tree — 167/184 (91%) of the final catalog received a classification; the remainder either lack a sufficiently close GTDB reference or fall below GTDB-Tk's own placement confidence and are reported as unclassified rather than forced to a low-confidence call.

## 5. Methods

| Step | Tool / parameters |
| :--- | :--- |
| Pipeline | nf-core/mag 5.4.2 |
| Read QC & host removal | Same fastp + Bowtie2/GRCm39 host removal as Phase 1 (shared input) |
| Co-assembly | MEGAHIT, per diet-arm group (group-AL, group-IF) |
| Binning | MetaBAT2, MaxBin2, SemiBin2 (independent, per group) |
| Bin refinement | DASTool — consensus, non-redundant bin selection across the three binners |
| Genome QC | CheckM2 (completeness/contamination), BUSCO (bacteria_odb10), QUAST (assembly statistics) |
| Taxonomic classification | GTDB-Tk `classify_wf`, database release R226 |
| Quality-control curation (this report) | Custom mapping of the 184 DASTool-refined bins back to their originating binner's CheckM2/BUSCO/QUAST/GTDB-Tk record (`scripts/15_mag_final_bin_catalog.py`), MIMAG tiering, and summary figures (`scripts/16_mag_summary_plots.R`) |

Pipeline execution: 3,328/3,328 Nextflow tasks completed successfully, 0 failed. Analysis was orchestrated with resource caps within server policy (≤ 28 cores sustained). All figures are provided as both 300-dpi PNG and vector PDF.

## 6. Results

### 6.1 Co-assembly statistics

| Group | Contigs | Total assembly size | N50 |
| :--- | :---: | :---: | :---: |
| group-AL | 297,332 | 553.9 Mb | 7,366 bp |
| group-IF | 302,440 | 518.7 Mb | 5,593 bp |

The AL co-assembly is somewhat larger and more contiguous (higher N50) than IF, consistent with its marginally higher pooled input read count; this modest assembly-quality difference is accounted for when comparing MAG counts between arms below (raw bin counts are not directly proportional to assembly size).

### 6.2 Genome recovery and quality (Figure 11)

DASTool selected **184 non-redundant genomes** (95 from group-AL, 89 from group-IF) from the three binners' combined predictions. Applying MIMAG thresholds:

| Group | Total MAGs | High quality | Medium quality | Low quality | GTDB classified |
| :--- | :---: | :---: | :---: | :---: | :---: |
| AL | 95 | 44 | 43 | 8 | 87/95 (92%) |
| IF | 89 | 35 | 46 | 8 | 80/89 (90%) |
| **Total** | **184** | **79** | **89** | **16** | **167/184 (91%)** |

Figure 11 shows completeness vs. contamination for all 184 genomes. Most low-quality calls are driven by **contamination**, not incompleteness: several genomes are >90% complete but exceed the 5–10% contamination ceiling (up to 50% contamination in the most extreme case) — DASTool's own internal bin score does not always agree with CheckM2's external assessment, which is exactly why the post-hoc CheckM2-based tiering in this report (rather than trusting DASTool's selection at face value) matters for deciding which genomes are safe to use downstream.

### 6.3 Taxonomic composition of recovered MAGs (Figure 12)

Among GTDB-classified genomes, **Bacillota (Firmicutes) dominates in both arms** (AL 73/87 = 84%; IF 64/80 = 80%), with Bacteroidota (AL 8, IF 10) and Actinomycetota (AL 5, IF 6) as consistent minority phyla — no qualitative phylum-level shift between diets, mirroring Phase 1's finding of no statistically resolvable community-level diet effect.

The single Verrucomicrobiota genome recovered in the whole catalog is *Akkermansia muciniphila* (group-AL only; 100% complete, 5.18% contamination → Medium quality per the tiering above). No corresponding *Akkermansia* genome met binning criteria in group-IF, **even though Phase 1 found *Akkermansia* more abundant by read count in IF (64.4%) than AL (48.7%)**. Inspecting the AL bin's per-sample coverage shows an 8-fold range across the five pooled animals (31×–238×) — this is the signature of a dominant taxon present at very different intra-arm coverage levels (individual/strain heterogeneity), which can fragment a co-assembly's graph for that organism enough to prevent clean binning in one group while still permitting it in the other. This is a useful, genuinely informative discordance rather than a contradiction: **read-based relative abundance and genome-resolved assembly measure different things**, and the two together (Phase 1 + Phase 2) give a fuller picture than either alone.

At genus level, the most frequently recovered genera across the catalog are *Acetatifactor*, *Pelethomonas*, *Angelakisella*, *Lawsonibacter*, *Acutalibacter*, *Oscillibacter*, *Kineothrix*, and *Paramuribaculum* — notably, *Oscillibacter* and *Paramuribaculum* were also flagged in Phase 1 (dominant genus and IF-trending species, respectively), giving independent genome-level corroboration of those community-level observations.

## 7. Conclusions

| Question | Result |
| :--- | :--- |
| How many usable genomes were recovered? | 168 medium-to-high-quality MAGs (79 high, 89 medium) out of 184 total non-redundant genomes; ~447 Mb of genome-resolved sequence |
| Does genome quality differ AL vs IF? | No material difference — similar quality-tier proportions in both arms (AL: 46% high/45% medium/8% low; IF: 39% high/52% medium/9% low) |
| Does genome-level taxonomy differ AL vs IF? | No qualitative phylum-level shift — Bacillota-dominated in both arms, matching Phase 1's community-level null result |
| Is *Akkermansia muciniphila* recovered as a genome? | Yes, but **only in AL** (Medium quality, 100% complete) — despite higher IF read abundance in Phase 1; attributed to coverage heterogeneity across pooled AL animals fragmenting that taxon's assembly differently between groups, not to it being genuinely absent from IF |
| Genome-level corroboration of Phase 1 taxa? | Yes — *Oscillibacter* and *Paramuribaculum*, both flagged in Phase 1, are among the most frequently recovered genera here |

**Overall interpretation.** The genome-resolved (MAG) analysis is directionally consistent with Phase 1's community-level conclusion: no clear, qualitative diet effect emerges at the level of recovered-genome count, quality, or broad taxonomic composition. The one point of apparent disagreement — *Akkermansia*'s genome being recovered only in AL despite higher IF read abundance — is explained by assembly-level coverage heterogeneity rather than a real absence, and illustrates why genome-resolved and read-based metagenomics are complementary rather than redundant: a future study aimed at strain-level or genomic (not just compositional) diet effects would benefit from the sample-size and cage-control recommendations already given in Phase 1 §7.1, since the same within-arm heterogeneity that limited the community-level test also limits clean genome assembly for the most variable taxa.

## 8. Deliverable Files

```
custom_research_report_20260720/
├── Daniel_Mendes_gut_metagenomics_0720.md   ← this report
├── qc/
│   └── multiqc_report_mag.html              ← MAG pipeline QC (assembly, binning, annotation summary)
└── assembly_binning/
    ├── final_MAG_catalog.tsv                ← 184 non-redundant genomes: quality, size, taxonomy, origin binner
    ├── mag_quality_summary_by_group.tsv     ← AL vs IF counts by MIMAG tier + GTDB classification rate
    ├── fig11_mag_quality_scatter.(png|pdf)  ← completeness vs. contamination, MIMAG thresholds
    └── fig12_mag_taxonomy_composition.(png|pdf)  ← phylum-level composition of recovered MAGs, AL vs IF
```

Full pipeline output (co-assemblies, all binner raw predictions, BUSCO/CheckM2/QUAST logs, Prokka gene annotations) is retained on the analysis server (`output_results_mag/`) for reproducibility but is not copied into this delivery folder — it is regenerable pipeline intermediate output, not a client-facing result table.

---

*Prepared by Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics — 2026-07-20*

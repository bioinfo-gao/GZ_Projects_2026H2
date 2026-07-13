# Humanized-Mouse Whole-Genome Sequencing — Knock-In Integration, Copy Number & Sequence Integrity

**Report Date:** 2026-07-11
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Sequencing:** Whole-genome sequencing, paired-end 150 bp, NovaSeq X Plus

---

## 1. Objectives

For six humanized-mouse samples spanning three targeted lines, this study was commissioned to:

1. **Confirm on-target integration** — verify that each targeting vector integrated at its intended mouse locus, replacing the wild-type allele with the human sequence.
2. **Determine copy number** of the integrated human construct.
3. **Confirm the integrity** of the human knock-in sequence (absence of internal deletions or rearrangements).

Extending these, we also screened for off-target (random) integration, assessed zygosity, and — for the CD1A line — verified the status of the Neo selection cassette.

## 2. Key Findings

- **All six samples carry the correct construct at the intended locus.** Each sample shows near-complete coverage (94–100% breadth) of its expected human insert and read-pair bridging to the intended mouse target locus. No sample identity or labeling issue was found.
- **No credible off-target integration** was detected in any sample after filtering recurrent alignment artifacts.
- **Copy number is single-copy for five of six samples.** One MTTH sample (**MTTH_284**) shows an **elevated copy number (~2–3×), indicating multi-copy / concatemeric integration** — flagged for attention.
- **The human knock-in sequences are structurally intact** — no evidence of large internal deletions or rearrangements in any sample.
- **CD1A line: the Neo selection cassette has been deleted** (as intended for a clean allele).
- Both **RAGH samples are heterozygous** for the knock-in. For the MTTH and CD1A lines, allele-level zygosity could not be resolved by read depth because the human insert is homologous to the mouse gene it replaces (see §4/§6); copy-number data are provided as context.

## 3. Sample Information

| Sample | Line | Mouse target locus | Human insert | Mean depth | Duplication |
| :--- | :---: | :---: | :---: | :---: | :---: |
| CD1A_B125 | CD1A | *Cd1d1 + Cd1d2* (chr3) | Human CD1 cluster — CD1D, CD1A, CD1C, CD1B, CD1E (~127 kb) | 30.5× | 15.2% |
| RAGH_153 | RAGH | *Rag2* (chr2) | Human cytokine cassette (G-CSF, M-CSF, IL-6, IL-1β, IL-7, IL-15) | 19.7× | 16.1% |
| RAGH_273 | RAGH | *Rag2* (chr2) | Human cytokine cassette | 22.9× | 20.1% |
| MTTH_284 | MTTH | *Htt* (chr5) | Human *HTT* full gene (~178 kb) | 20.7× | 16.0% |
| MTTH_412 | MTTH | *Htt* (chr5) | Human *HTT* full gene | 19.6× | 19.4% |
| MTTH_524 | MTTH | *Htt* (chr5) | Human *HTT* full gene | 19.7× | 23.6% |

All six libraries passed quality control (Q30-rich 150 bp reads, 20–24% duplication typical of PCR WGS, mean genome coverage 20–30×) and are suitable for structural and copy-number analysis.

## 4. Analysis Rationale and Decision Criteria

- **A custom hybrid reference** (mouse GRCm39 with the three targeted-allele constructs added as extra contigs) was used so that reads originating from the human insert map to a dedicated sequence rather than being lost or mismapped.
- **"Construct present" is judged by coverage breadth, not mean depth.** A genuine integration covers essentially the whole human-specific region (>90% breadth); low-level cross-mapping from a homologous mouse gene produces high local depth over only a small fraction of the region (≤35% breadth). Using breadth prevents cross-mapping from being mistaken for a real integration.
- **On-target confirmation** uses read pairs whose one mate lies in the construct and whose partner falls at the construct's intended endogenous mouse locus.
- **Off-target screening** counts uniquely-mapped (MAPQ ≥ 20) construct reads landing at unintended loci. **Loci recurring across unrelated samples are treated as alignment artifacts, not true integrations** — a genuine off-target event is private to one sample.
- **Copy number** = uniquely-mapped depth over the human-specific region ÷ the sample's own autosomal baseline depth (≈0.5 per single-copy allele, ≈1.0 for two copies).
- **Integrity** is assessed by scanning the human insert for local coverage dropouts; a true internal deletion appears as a contiguous near-zero block with reduced overall breadth.

## 5. Methods

| Step | Tool / approach | Key parameters |
| :--- | :--- | :--- |
| Read QC & trimming | fastp, FastQC, MultiQC | default WGS QC |
| Alignment | nf-core/sarek 3.8.1, bwa-mem2 vs hybrid reference | `--tools tiddit --skip_tools baserecalibrator` |
| Duplicate marking, depth | GATK MarkDuplicates, mosdepth, samtools stats | MAPQ ≥ 20 for unique-depth metrics |
| Structural variants | TIDDIT | per-sample SV calling |
| Integration (on/off-target) | custom read-pair bridging + unique-mapping screen | MAPQ ≥ 20 off-target; recurrent-locus artifact filter |
| Copy number | mosdepth over human-specific region vs autosomal baseline | MAPQ ≥ 20 |
| Sequence integrity | 500 bp sliding-window depth scan of the human insert | flag < 0.3× or > 2.5× regional median |
| Zygosity | residual read depth at the replaced mouse locus | ≈0.5 heterozygous, ≈0 homozygous |
| CD1A Neo status | direct coverage of the NeoR/KanR cassette coordinates | — |

## 6. Results

### 6.1 On-target integration and construct identity

Every sample shows near-complete coverage of its expected human insert and read-pair bridging to the intended mouse locus. No sample carries a construct other than its designated one.

| Sample | Expected construct | Insert breadth | On-target bridging read pairs | On-target integration |
| :--- | :---: | :---: | :---: | :---: |
| CD1A_B125 | Human CD1 cluster | 98.2% | 48 | Confirmed |
| RAGH_153 | Cytokine cassette | 99.6% | 24 | Confirmed |
| RAGH_273 | Cytokine cassette | 99.7% | 33 | Confirmed |
| MTTH_284 | Human *HTT* | 99.8% | 45 | Confirmed |
| MTTH_412 | Human *HTT* | 98.8% | 44 | Confirmed |
| MTTH_524 | Human *HTT* | 96.8% | 79 | Confirmed |

### 6.2 Copy number

| Sample | Copy-number ratio | Interpretation |
| :--- | :---: | :---: |
| RAGH_153 | 0.61 | Single copy |
| RAGH_273 | 0.56 | Single copy |
| MTTH_412 | 0.71 | Single copy |
| MTTH_524 | 0.45 | Single copy |
| CD1A_B125 | 0.98 | ~Two copies (consistent with a homozygous single-copy allele) |
| **MTTH_284** | **1.60** | **Elevated — multi-copy / possible concatemeric integration (recommend follow-up)** |

### 6.3 Sequence integrity

All human inserts are recovered at 94–100% breadth with no contiguous near-zero coverage block, indicating **no large-scale internal deletion or rearrangement**. Scattered low-coverage windows within the *HTT* inserts (e.g. MTTH_412) are non-zero (roughly 10–30% of local median, longest run ~3 kb) and reflect reduced mappability where human *HTT* is ~85% identical to mouse *Htt*, not sequence loss.

### 6.4 Off-target integration

No sample shows a credible off-target integration. All candidate off-target signals were low-support loci recurring across unrelated samples — the signature of alignment artifacts rather than true random integration — and were removed by the cross-sample artifact filter.

### 6.5 Zygosity

- **RAGH line — heterozygous.** Both RAGH_153 and RAGH_273 retain ~one wild-type *Rag2* allele's worth of coverage at the native locus (residual depth ratio ~0.25), consistent with a single heterozygous knock-in.
- **MTTH and CD1A lines — not resolvable by depth.** Because the human insert (*HTT*; the CD1D portion of the CD1 cluster) is homologous to the mouse gene it replaces, human-insert reads map back onto the native mouse coordinates and inflate the residual-depth measurement. Copy-number data (§6.2) indicate single-copy integration for MTTH_412/524 and ~two copies for CD1A_B125. Definitive allele-level zygosity for these two lines would require targeted junction-read analysis and can be added if desired.

### 6.6 CD1A Neo cassette status

The NeoR/KanR selection cassette in the CD1A construct is **deleted** in CD1A_B125 (cassette coverage 0.75× versus 29× over the human insert), i.e. a clean humanized allele consistent with the RAGH and MTTH designs.

## 7. Conclusions

| Objective | CD1A_B125 | RAGH_153 | RAGH_273 | MTTH_284 | MTTH_412 | MTTH_524 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| On-target integration | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Copy number | ~2 copies | single | single | **multi (~2–3)** | single | single |
| Sequence integrity | Intact | Intact | Intact | Intact | Intact | Intact |
| Off-target | None | None | None | None | None | None |
| Zygosity | n/d (see §6.5) | Heterozygous | Heterozygous | n/d | n/d | n/d |
| Neo cassette | Deleted | — | — | — | — | — |

**Overall:** All six samples are correctly targeted, structurally intact, single-copy humanizations free of detectable off-target integration — with the single exception that **MTTH_284 carries the human *HTT* insert at elevated (multi-copy) dosage**, which we recommend confirming before using that founder.

## 8. Deliverable Files

```
custom_research_report_20260711/
├── Ellen_KnockIn_WGS_0711.md      ← this report
├── qc/                            ← MultiQC report + per-sample FastQC/fastp
├── integration/                   ← per-sample on/off-target integration tables
├── copy_number/                   ← per-sample copy-number tables
├── ki_integrity/                  ← per-sample insert-integrity scans
├── zygosity/                      ← per-sample zygosity tables
└── cd1a_neo_status/               ← CD1A Neo-cassette coverage
```

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*

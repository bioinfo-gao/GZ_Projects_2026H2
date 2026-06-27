# Mouse RNA-seq Project — Candidate Bacterial/Mycoplasma Contamination Sources

**Report Date**: 2026-06-27
**Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform**: Linux HPC server

---

## 1. Purpose of This Note

The client's follow-up question specifically asked for a **breakdown of bacterial/Mycoplasma reads** behind the "Bacterial 16S"-dominant rRNA signal seen in mouse_28, mouse_29, and mouse_46 (main report, Section 6.2/7.2). The main report (`custom_research_report_20260627/Mouse_rRNA_QC_0627.md`) answered this conservatively using only the SortMeRNA rRNA-database data, noting that the generic 16S match cannot identify species.

This note goes one step further and reports what a **species-level screen** of the same data actually found — framed explicitly as **candidate/possible contamination sources**, not confirmed identifications, for the reasons explained in Section 4 below.

---

## 2. Method

Reads that failed to align to the mouse genome (GRCm39) in two representative samples were screened against Kraken2's Standard-8 database (RefSeq archaea + bacteria + viral + plasmid + human + UniVec_Core):

- **mouse_29** — representative of the severe-dimer / Bacterial-16S-dominant group (mouse_28 shares the same library profile and was not screened separately to avoid redundant testing).
- **mouse_46** — the distinct, mixed-pattern outlier with near-full-length reads but the lowest mapping rate in the batch.

For each candidate organism, the table below reports the **higher of the two samples' values** — i.e., the upper bound we observed, not an average — to keep the framing conservative and avoid understating the signal.

---

## 3. Candidate Contamination Sources Observed

| Candidate source | Reads observed (higher of the two samples) | % of that sample's 10M-pair QC subsample | What it typically indicates |
| :--- | :---: | :---: | :---: |
| *Mesomycoplasma hyorhinis* (Mycoplasma) | 3,898 (mouse_46) | 0.039% | Common, hard-to-detect cell-culture/tissue contaminant |
| *Vibrio furnissii* | 12,019 (mouse_46) | 0.120% | Environmental/water-associated bacterium; commonly a reagent or lab-water background signal, not a biological contaminant of the tissue itself |
| *Legionella pneumophila* | 391 (mouse_46) | 0.004% | Water-system-associated bacterium; same likely source as above |
| *Bacillus amyloliquefaciens* | 2,134 (mouse_46) | 0.021% | Common soil/environmental bacterium, frequent low-level lab background |
| *Actinoalloteichus* sp. | 2,623 (mouse_46) | 0.026% | Environmental actinobacterium, frequent low-level lab background |
| All bacteria combined | 94,084 (mouse_46) | 0.94% | Sum of all bacterial hits, all organisms combined |

**For scale:** even the single largest specific hit — *Mycoplasma* — amounts to about **4 reads in every 10,000 sequenced**, in the more affected of the two samples tested. The other named organisms are smaller still.

---

## 4. Why These Are "Possible" Sources, Not Confirmed Contamination

1. **Low-abundance species calls from short RNA-seq reads carry inherent uncertainty.** A 150 bp read matching a reference at this depth is suggestive, not definitive — species-level calls on this kind of data are best treated as leads to investigate, not laboratory-confirmed findings.
2. **The reference database used does not include the mouse genome.** As a result, a large share of these same reads — 11.6% (mouse_46) to 12.4% (mouse_29) of all reads screened — were assigned to *Homo sapiens*, simply because human is the closest match available in this database for sequence shared between mouse and human (a normal consequence of evolutionary conservation, not necessarily literal human DNA). This is the same caveat that applies, in principle, to the smaller bacterial calls: without a mouse reference in the comparison set, some "bacterial-looking" reads could in fact be unrecognized mouse sequence with a coincidental partial match.
3. **Only 2 of the 9 samples were screened this way**, and only their already-unaligned reads (a small subset of each library) — this is a spot-check, not a systematic survey of the batch.

For these reasons, the table in Section 3 should be read as: *"these are the specific organisms our screen pointed to, at the stated low frequency, in the worse of the two samples tested" — not "we have confirmed the sample is contaminated with these organisms."*

---

## 5. Conclusion and Recommendation

- The species-level screen is **consistent with** the main report's conclusion (Section 7.2) that the 15.6–37.7% "Bacterial 16S" figure from the rRNA tool is an overestimate: actual candidate bacterial signal, even taken at face value, tops out under 1% of reads — roughly 30–60× smaller than the rRNA tool's number.
- None of the candidate sources found are large enough to explain the batch's rRNA contamination or mapping-rate problems; those remain attributable to the adapter-dimer and rRNA-depletion issues described in the main report.
- The one finding worth a precautionary follow-up, independent of its small size, is ***Mycoplasma*** — it is a well-known contaminant of cultured cells/tissue that is invisible by eye or under a standard microscope. We recommend a routine, low-cost PCR-based Mycoplasma test on the original source material as good practice, not because this data proves contamination, but because the candidate signal — however small — is specific enough to be worth ruling out at the bench.

---

*Report prepared by:*

**Zhen Gao, PhD**
Principal Bioinformatics Scientist
Athenomics

# sgRNA Guide Reference — Study A (all 9 guides)

Product: CRISPRevolution sgRNA EZ Kit (Modified), 1.5 nmol each. Sequences were supplied by the client as
images (RNA); converted here to DNA spacer sequences. Cut sites are located **by sequence** against
GRCm39 (GENCODE vM35), not by the coordinates embedded in the product names — see the note below.

## Guides as supplied

| Target | Guide ID / name | Spacer (DNA, 5'→3') | Used in |
| :--- | :--- | :--- | :---: |
| Pten | Pten-32799878 | `GGTGGGTTATGGTCTTCAAA` | B1TP + B2TP |
| Pten | Pten-32799895 | `TGATAAGTTCTAGCTGTGGT` | B1TP + B2TP |
| Pten | Pten-32799899 | `GGTTTGATAAGTTCTAGCTG` | B1TP + B2TP |
| Brca1 | (guide 1) | `GGTTCCGGTAGCCCACGCTC` | B1TP |
| Brca1 | (guide 2) | `GGCGTCGATCATCCAGAGCG` | B1TP |
| Brca1 | (guide 3) | `TTCTTGTGAGCGTTTGAATG` | B1TP |
| Brca2 | Brca2+150529497 | `GATAAGCCTCAATTGGTTTG` | B2TP |
| Brca2 | Brca2−150529492 | `AAAGCTCCTCAAACCAATTG` | B2TP |
| Brca2 | Brca2−150529524 | `AGGTTCAGAATTGTATGGGG` | B2TP |

## Brca2 guide verification (performed 2026-07-16, before any downstream use)

Each Brca2 spacer maps to **exactly one** site in GRCm39, on the strand indicated by the client's own
naming, and each carries a canonical NGG PAM. This confirms the guides are genuine and that the
transcription from the supplied image is correct.

| Guide | GRCm39 location (1-based) | Strand | Protospacer + PAM | Predicted cut site |
| :--- | :---: | :---: | :---: | :---: |
| Brca2+150529497 | chr5:150,452,945–150,452,964 | + | `GATAAGCCTCAATTGGTTTG` + `AGG` | chr5:150,452,961 |
| Brca2−150529492 | chr5:150,452,954–150,452,973 | − | `AAAGCTCCTCAAACCAATTG` + `AGG` | chr5:150,452,957 |
| Brca2−150529524 | chr5:150,452,986–150,453,005 | − | `AGGTTCAGAATTGTATGGGG` + `GGG` | chr5:150,452,989 |

All three cut sites lie within **Brca2 exon 3 (CDS)** of the Ensembl-canonical / CCDS transcript
`ENSMUST00000044620.11` (*Brca2-201*, CCDS39411.1).

**Design:** guides 1 and 2 overlap and cut the same position from opposite strands (4 bp apart);
guide 3 cuts ~30 bp downstream. The three cuts fall within a **33 bp window**
(150,452,957–150,452,989), a multi-guide design predicting either small indels or excision of the
intervening fragment. The observed 31 bp deletion in RO_B2TP is exactly that excision.

## Note on the coordinates in the guide IDs

The numbers embedded in the Brca2 product names (e.g. `150529497`) **are not GRCm39 coordinates**:
they lie ~76.5 kb away from where the guides actually map on GRCm39, so using them directly against
GRCm39 retrieves the wrong locus. The offset is most consistent with GRCm38/mm10 numbering, although
this was not verified directly (no GRCm38 reference was used in this project). Because every guide
here was located by sequence, no result is affected.

---

*Zhen Gao, PhD — Principal Bioinformatics Scientist, Athenomics*

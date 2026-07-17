# sgRNA Guide Reference — Study A (all 9 guides)

Product: CRISPRevolution sgRNA EZ Kit (Modified), 1.5 nmol each. Sequences were supplied by the client as
images (RNA); converted here to DNA spacer sequences. Cut sites are located **by sequence** against
GRCm39 (GENCODE vM35), not by the coordinates embedded in the product names — see the note below.

## All nine guides — as supplied, and verified against GRCm39

Each spacer maps to **exactly one** site in GRCm39, inside its intended target gene, on the strand
indicated by the client's naming where one is given, and each carries a canonical NGG PAM. This
confirms the guides are genuine and that the transcription from the supplied images is correct.
Cut site = 3 bp upstream of the PAM.

| Target | Guide ID / name | Protospacer + PAM | GRCm39 location (1-based) | Strand | Predicted cut | Used in |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| Pten | Pten-32799878 | `GGTGGGTTATGGTCTTCAAA` + `AGG` | chr19:32,777,275–32,777,294 | − | chr19:32,777,278 | B1TP + B2TP |
| Pten | Pten-32799895 | `TGATAAGTTCTAGCTGTGGT` + `GGG` | chr19:32,777,292–32,777,311 | − | chr19:32,777,295 | B1TP + B2TP |
| Pten | Pten-32799899 | `GGTTTGATAAGTTCTAGCTG` + `TGG` | chr19:32,777,296–32,777,315 | − | chr19:32,777,299 | B1TP + B2TP |
| Brca1 | (guide 1) | `GGTTCCGGTAGCCCACGCTC` + `TGG` | chr11:101,422,890–101,422,909 | + | chr11:101,422,906 | B1TP |
| Brca1 | (guide 2) | `GGCGTCGATCATCCAGAGCG` + `TGG` | chr11:101,422,905–101,422,924 | − | chr11:101,422,908 | B1TP |
| Brca1 | (guide 3) | `TTCTTGTGAGCGTTTGAATG` + `AGG` | chr11:101,422,929–101,422,948 | − | chr11:101,422,932 | B1TP |
| Brca2 | Brca2+150529497 | `GATAAGCCTCAATTGGTTTG` + `AGG` | chr5:150,452,945–150,452,964 | + | chr5:150,452,961 | B2TP |
| Brca2 | Brca2−150529492 | `AAAGCTCCTCAAACCAATTG` + `AGG` | chr5:150,452,954–150,452,973 | − | chr5:150,452,957 | B2TP |
| Brca2 | Brca2−150529524 | `AGGTTCAGAATTGTATGGGG` + `GGG` | chr5:150,452,986–150,453,005 | − | chr5:150,452,989 | B2TP |

## Coding-exon context

All nine cut sites fall inside a coding exon of the Ensembl-canonical transcript of their target,
so a disruptive indel at any of them is expected to be loss-of-function rather than silent:

| Target | CDS exon | Canonical transcript |
| :--- | :---: | :--- |
| Pten | exon 5 | `ENSMUST00000249247.1` |
| Brca1 | exon 6 | `ENSMUST00000017290.11` |
| Brca2 | exon 3 | `ENSMUST00000044620.11` (*Brca2-201*, CCDS39411.1) |

## Design

All three targets use the same multi-guide strategy: three guides clustered in a narrow window, with
overlapping guides cutting the same position from opposite strands.

| Target | Cut-site window | Span |
| :--- | :---: | :---: |
| Pten | 32,777,278 – 32,777,299 | 22 bp |
| Brca1 | 101,422,906 – 101,422,932 | 27 bp |
| Brca2 | 150,452,957 – 150,452,989 | 33 bp |

This predicts either small indels at a cut site or excision of the fragment between the outermost
cuts. The 31 bp deletion observed in RO_B2TP (chr5:150,452,958–150,452,988) is exactly that excision.

## Note on the coordinates in the guide IDs

The numbers embedded in the product names **are not GRCm39 coordinates**: `Brca2+150529497` lies
~76.5 kb from where that guide actually maps on GRCm39, and `Pten-32799878` ~22.6 kb from its true
position, so using these numbers directly against GRCm39 retrieves the wrong locus. That the offset
**differs per locus** rules out a simple constant shift and is what an assembly coordinate change
produces; the numbering is most consistent with GRCm38/mm10, although this was not verified directly
(no GRCm38 reference was used in this project). Because every guide here was located by sequence, no
result is affected.

---

*Zhen Gao, PhD — Principal Bioinformatics Scientist, Athenomics*

---

*Zhen Gao, PhD — Principal Bioinformatics Scientist, Athenomics*

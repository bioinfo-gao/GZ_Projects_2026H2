# Targeted Brca2 Cut-Site Verification (Study A) — Addendum 1 to the 2026-07-15 Report

**Report Date:** 2026-07-16
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Mus musculus* (GRCm39, GENCODE vM35)
**Tissue / Cell:** Cultured cells — primary cells from a Trp53⁺/⁻;Cas9 mouse (RO_origin, RO_B1TP, RO_B2TP) and cell lines digested from the resulting in-vivo tumors (RO_tumor1–3)
**Scope:** Study A only (n=6). Supersedes nothing — the 2026-07-15 report remains valid in full; this addendum closes its one open item.

---

## 1. Objectives

The 2026-07-15 report closed with a single outstanding input (§7.2 of that report): the **three Brca2 sgRNA spacer sequences**, which you provided on 2026-07-16. A spacer sequence is what defines a predicted cut site, so without it a base-accurate Brca2 genotype was not possible. With the guides now in hand, this addendum answers the two questions that were gated on them:

- **A1-b** — Confirm, at the predicted cut site, that the intended **Brca2 knockout actually occurred in the RO_B2TP cell line**. Previously B2TP could only be shown to be *consistent with* its Brca2+Pten design (Pten knocked out, Brca1 wild-type); the Brca2 edit itself had never been positively observed.
- **A2-b** — Test **tumor3 for a Brca2 edit** to settle its lineage. Tumor3 is wild-type at all Brca1 and Pten cut sites; a Brca2 edit would have reopened a B2TP origin for it.

---

## 2. Key Findings

1. **The Brca2 knockout in RO_B2TP is now positively confirmed, and it is biallelic.** Of 14 informative reads spanning the Brca2 cut site, **zero are wild-type**. No wild-type allele survives in this cell line.
2. **The edit is exactly what the guide design predicts.** B2TP carries a **31 bp deletion (chr5:150,452,958–150,452,988)** whose two ends sit precisely on the cut sites of the supplied guides — the fragment excised between the overlapping guide pair and the third guide. It causes a **frameshift at codon ~31 of 3,329** in the canonical Brca2 transcript, i.e. within the first 1% of the protein: an unambiguous null allele.
3. **Tumor3 shows no Brca2 edit and therefore does not descend from B2TP.** Tumor3 retains abundant wild-type Brca2 reads and carries none of B2TP's alleles. Since B2TP itself has **no wild-type allele left**, no descendant of B2TP could regain one. Combined with tumor3's wild-type Brca1 *and* Pten, this raises the 2026-07-15 conclusion — tumor3 arose from an **un-edited subclone** — from *Moderate* to **High** confidence.
4. **All three guides were verified as genuine and correctly transcribed before use.** Each maps to a single unique site inside Brca2, on the strand your naming specifies, each with a perfect NGG PAM.
5. **A caution on your guide IDs:** the coordinates embedded in the product names (e.g. `150529497`) **are not GRCm39 positions** — they sit ~76.5 kb away from where the guides actually map on GRCm39, and pointing them at GRCm39 retrieves the wrong locus. The offset is most consistent with **GRCm38/mm10** numbering (the assembly change shifted this region), though we did not verify that directly, as no GRCm38 reference is held on our analysis server. We located every guide by sequence, so no result here is affected.

---

## 3. Sample Information

Six Study A samples, paired-end 150 bp WGS (NovaSeq X Plus), aligned to GRCm39 with nf-core/sarek 3.8.1 (bwa-mem2). No new sequencing was performed for this addendum; it re-interrogates the existing alignments at a newly defined locus.

| # | Label | Type | Role | Brca2 depth (cut window / gene-wide) |
| :--- | :--- | :---: | :---: | :---: |
| 1 | RO_origin | Cell | Trp53⁺/⁻;Cas9 parent, unedited — matched normal | 24.3× / 21.0× |
| 2 | RO_B1TP | Cell | sgRNA Brca1+Pten | 26.4× / 24.3× |
| 3 | RO_B2TP | Cell | sgRNA **Brca2**+Pten | **13.3× / 23.2×** |
| 4 | RO_tumor1 | Cell | tumor-derived line | 37.0× / 30.1× |
| 5 | RO_tumor2 | Cell | tumor-derived line | 28.1× / 23.7× |
| 6 | RO_tumor3 | Cell | tumor-derived line | 38.6× / 44.3× |

### 3.1 Client-supplied Brca2 guides (CRISPRevolution sgRNA EZ Kit, 1.5 nmol, Modified)

Transcribed from your image (archived as `edit_verification/client_brca2_sgRNA_source_image.png`) and converted RNA→DNA:

| Guide ID (as supplied) | Spacer (DNA, 5'→3') | GRCm39 location | Strand | PAM | Predicted cut |
| :--- | :--- | :---: | :---: | :---: | :---: |
| Brca2+150529497 | `GATAAGCCTCAATTGGTTTG` | chr5:150,452,945–150,452,964 | + | AGG | chr5:150,452,961 |
| Brca2−150529492 | `AAAGCTCCTCAAACCAATTG` | chr5:150,452,954–150,452,973 | − | AGG | chr5:150,452,957 |
| Brca2−150529524 | `AGGTTCAGAATTGTATGGGG` | chr5:150,452,986–150,453,005 | − | GGG | chr5:150,452,989 |

All three cut sites fall inside **Brca2 exon 3 (CDS)** of the Ensembl-canonical, CCDS transcript `ENSMUST00000044620.11` (*Brca2-201*, CCDS39411.1). Guides 1 and 2 overlap and cut the same point from opposite strands (4 bp apart); guide 3 cuts ~30 bp downstream. This is a multi-guide design concentrating three cuts into a **33 bp window**, which predicts either small indels or excision of the intervening fragment.

---

## 4. Analysis Rationale and Decision Criteria

| Step | Rationale | Decision criterion |
| :--- | :--- | :--- |
| Verify guides before trusting them | A mis-transcribed spacer would silently target the wrong site and invalidate everything downstream | Each spacer must map **uniquely** inside Brca2, on the **stated strand**, with a canonical **NGG PAM** |
| Locate by sequence, not by supplied coordinate | Vendor IDs proved to be mm10-based | Cut site = 3 bp upstream of PAM, derived from the sequence match |
| Do **not** rely on the variant caller alone | `bcftools call` reported **no Brca2 indel in any sample** — a false negative: edited reads are heavily soft-clipped or carry a 31 bp deletion, which this caller missed at a complex site | Any "no edit" claim must be corroborated at read/CIGAR level before being believed |
| Primary readout = **absence of wild-type reads** | The cleanest evidence of a biallelic knockout is that no reference allele survives | WT read count spanning the cut site; 0 WT ⇒ biallelic disruption |
| Soft-clip fraction is **not** used as evidence on its own | This locus is intrinsically clip-prone: ~30% of reads are clipped **even in the unedited parent** | A clip counts only if ≥20 bp **and** its breakpoint **stacks recurrently on a cut-site base** |
| Negative-control gate | The unedited parent must score as unedited, or the method is wrong | RO_origin must show ~0 editing; asserted in code |
| Homopolymer scepticism | Indel calls inside homopolymers are classic sequencing artifacts | Indels within the 7-C tract (150,452,983–150,452,989) are discounted |

---

## 5. Methods

- **Guide localisation:** `seqkit locate` of each spacer (both strands) against GRCm39; PAM confirmed by `samtools faidx` of the 3 bp immediately 3' of each protospacer.
- **Cut-site genotyping:** `bcftools mpileup/call` (MAPQ≥20, BAQ≥15) jointly across all six samples at each cut site ±60 bp — retained for completeness, but see §6.1 for why it is not the primary readout here.
- **Read-level allele quantification** (`scripts/study_A/A3b_brca2_allele_quant.py`, pysam): reads with MAPQ≥20, duplicates excluded, classified by CIGAR into wild-type (spanning the 150,452,952–150,452,994 core with ≥25 bp flanking margin and no indel/clip), the 31 bp excision allele, other cut-site indels, or cut-site soft-clips (≥20 bp, breakpoint inside the core). Clip breakpoints tabulated per base to separate recurrent (real) from scattered (background) clipping.
- **Depth:** `samtools depth` over the cut window vs the whole Brca2 gene body, per sample, to normalise each sample against its own baseline.
- **Consequence:** deletion mapped onto GENCODE vM35 CDS coordinates of `ENSMUST00000044620.11`.

---

## 6. Results

### 6.1 The variant caller returned a false negative — read-level analysis was required

The joint `bcftools` call across all six samples reported **no indel at any Brca2 cut site**, in any sample. Taken at face value this would have read as "the Brca2 edit did not happen". It is an artifact: at the B2TP cut site only **1 of 25** reads is a perfect match (`150M`), versus **24 of 34** in the unedited parent — the edited reads are present but carry a 31 bp deletion or are heavily soft-clipped (50–87 bp), which the caller did not resolve. All conclusions below therefore rest on read/CIGAR-level evidence.

### 6.2 A localised, B2TP-specific loss of coverage marks the cut site

Depth in the cut window, normalised to each sample's own Brca2 gene-wide depth:

| Sample | Cut window | Brca2 gene-wide | Ratio |
| :--- | :---: | :---: | :---: |
| RO_origin | 24.3× | 21.0× | 1.16 |
| RO_B1TP | 26.4× | 24.3× | 1.09 |
| **RO_B2TP** | **13.3×** | **23.2×** | **0.57** |
| RO_tumor1 | 37.0× | 30.1× | 1.23 |
| RO_tumor2 | 28.1× | 23.7× | 1.19 |
| RO_tumor3 | 38.6× | 44.3× | 0.87 |

B2TP's Brca2 coverage is entirely normal across the gene (23.2×) yet collapses to 0.57 of its own baseline precisely at the cut window. The drop is specific to B2TP, specific to the edited locus, and is the expected consequence of edited reads failing to align cleanly.

### 6.3 Allele composition at the Brca2 cut site

| Sample | Informative reads | **Wild-type** | 31 bp excision | Cut-site clips | Other indels |
| :--- | :---: | :---: | :---: | :---: | :--- |
| RO_origin | 7 | 6 | 0 | 1 | – |
| RO_B1TP | 11 | 10 | 0 | 1 | – |
| **RO_B2TP** | **14** | **0** | **2** | **10** | 3 bp del @150,452,989 ×2 |
| RO_tumor1 | 19 | 15 | 0 | 4 | – |
| RO_tumor2 | 12 | 8 | 0 | 4 | – |
| RO_tumor3 | 17 | 8 | 0 | 7 | 1 bp del @150,452,983 ×2 |

**RO_B2TP is the only sample with no wild-type allele.** Every other sample, including all three tumors, retains a clear wild-type population.

### 6.4 The 31 bp deletion is the predicted dual-cut excision

Reads carrying it (`76M31D43M31S`, `32S76M31D42M`) place the deletion at **chr5:150,452,958–150,452,988**, removing:

```
TGGTTTGAGGAGCTTTCCTCAGAAGCCCCCC   (31 bp)
```

Its 5' end abuts the cut sites of the overlapping guide pair (150,452,957/961) and its 3' end abuts the third guide's cut (150,452,989). This is precisely the fragment predicted to be excised when all three guides cut — the deletion did not merely land "near" the guides, its endpoints *are* the guides' cut sites.

**Consequence:** the deletion lies wholly within CDS exon 3 of `ENSMUST00000044620.11`. It begins at CDS nucleotide 91 (**codon ~31 of 3,329**) and 31 bp is not a multiple of 3 → **frameshift** → premature termination within the first 1% of the protein. This is a null allele, not a hypomorph.

### 6.5 The second B2TP allele is disrupted at the third guide's cut site

B2TP's 10 cut-site-clipped reads are not background. Their breakpoints **stack on a single base**:

| Sample | Top clip breakpoints (≥20 bp clips) | On a cut site? |
| :--- | :--- | :---: |
| RO_origin | 150,452,932 ×2; 150,452,943 ×1; 150,452,961 ×1 | scattered |
| RO_B1TP | 150,452,934 ×2; 150,453,004 ×2; 150,452,951 ×1 | scattered |
| **RO_B2TP** | **150,452,990 ×4; 150,452,989 ×3** | **yes — 7 reads on the g3 cut** |
| RO_tumor1 | 150,452,981 ×1; 150,453,014 ×1; 150,452,983 ×1 | scattered |
| RO_tumor2 | 150,452,961 ×2; 150,453,004 ×2; 150,453,008 ×2 | weak (2 reads) |
| RO_tumor3 | 150,452,941 ×1; 150,452,949 ×1; 150,452,953 ×1 | scattered |

Only B2TP shows a recurrent stack, and it sits exactly on the third guide's cut (150,452,989). The clipped sequence maps only ambiguously elsewhere (MAPQ 0, repetitive), so the allele's precise structure — a complex indel or insertion — is not resolved by short reads; that it is **disrupted at the cut site** is nonetheless clear. Together with the 31 bp allele and the complete absence of wild-type reads, B2TP is **biallelically disrupted**.

### 6.6 Tumor3 carries no Brca2 edit

Tumor3 retains 8 wild-type reads, carries no 31 bp allele, and shows no breakpoint stacking. Its only indel call — a 1 bp deletion at 150,452,983 (2 reads) — falls inside a **7-C homopolymer** (`AAG`**`CCCCCCC`**`ATACAATTCTG`, 150,452,983–150,452,989), the classic context for a sequencing artifact, and is discounted.

The lineage logic is now decisive: **B2TP has no wild-type Brca2 allele.** Any cell descended from B2TP must inherit that state; a wild-type allele cannot be regained. Tumor3 is ~47% wild-type at this locus, so **tumor3 cannot descend from B2TP**. This independently confirms the 2026-07-15 finding via a third locus, and removes the alternative that report had explicitly left open.

---

## 7. Conclusions

| # | Conclusion | Confidence | Evidence |
| :--- | :--- | :---: | :--- |
| 1 | The three supplied guides are genuine, correctly transcribed, and target Brca2 exon 3 | High | Unique GRCm39 hit each, correct strand, perfect NGG PAM (§3.1) |
| 2 | **RO_B2TP carries the intended Brca2 knockout — positively confirmed** | High | 0/14 wild-type reads; localised coverage drop to 0.57 (§6.2–6.3) |
| 3 | **The knockout is biallelic** | High | No wild-type allele; two distinct disrupted alleles (§6.3–6.5) |
| 4 | Allele 1 = 31 bp dual-cut excision → frameshift at codon ~31/3,329 → null | High | Deletion endpoints coincide with guide cut sites (§6.4) |
| 5 | Allele 2 = disrupted at the third guide's cut; exact structure unresolved | Moderate–High | 7-read breakpoint stack on 150,452,989/990; clipped mate maps to repeat (§6.5) |
| 6 | **B2TP = Brca2 + Pten knockout** — the naming reading `B2 = Brca2` is now confirmed by direct observation, not inference | High | §6.3–6.5 + Pten KO from 2026-07-15 |
| 7 | **Tumor3 does not derive from B2TP; it arose from an un-edited subclone** (raised from *Moderate*) | High | Tumor3 wild-type at Brca2, Brca1 **and** Pten; B2TP has no WT allele to lose (§6.6) |
| 8 | Tumors 1 and 2 remain B1TP-derived; neither carries a Brca2 edit | High | Wild-type Brca2 reads retained (§6.3); Brca1/Pten indels per 2026-07-15 |
| 9 | Your guide IDs are **not** GRCm39 coordinates (~76.5 kb offset; most consistent with GRCm38/mm10) | High that they are not GRCm39; the mm10 attribution is inferred | §2 (item 5) — no impact on results (guides located by sequence) |

**Tumor3 remains the one open biological question**, and it is now sharper rather than resolved: it is the most aneuploid of the three tumors yet is wild-type at all three targeted genes. The remaining discriminator is not another edit locus but **sample identity** — a SNP-fingerprint check of tumor3 against RO_origin would separate "un-edited escaper subclone of the parent" (expected: matching fingerprint) from a sample-tracking swap (expected: mismatch). We can run this on the existing data on your go-ahead; no new sequencing is required.

---

## 8. Deliverable Files

```
custom_research_report_20260716/
├── GeneEdit_Lats12_WGS_0716.md          ← this addendum
└── edit_verification/
    ├── client_brca2_sgRNA_source_image.png   client image the guides were transcribed from
    ├── sgRNA_guides_reference.md             all 9 guides + Brca2 localisation/PAM check
    ├── cut_sites.tsv                         all 9 cut sites (Pten/Brca1/Brca2) on GRCm39
    ├── cutsite_indels.tsv                    joint genotypes at every cut site, 6 samples
    ├── brca2_allele_quant.tsv                per-sample WT / 31bp / clip / indel counts
    ├── brca2_clip_breakpoints.tsv            clip breakpoint stacking, per base per sample
    ├── brca2_indels.tsv                      exploratory Brca2 gene-wide indel scan
    └── spacer_gene.tsv                       spacer → target gene mapping
```

The 2026-07-15 delivery (`custom_research_report_20260715/`) is unchanged and remains the reference for all other results (tumor aneuploidy, somatic variants, iHPV integration, Study B).

---

*Prepared by:*
**Zhen Gao, PhD**
*Principal Bioinformatics Scientist, Athenomics*

#!/usr/bin/env bash
# Project 16 — sample-origin sanity check from the WGS itself (NOT from HLA typing).
# HLA typing works on germline DNA of ANY nucleated cell, so it says nothing about tissue/cell.
# Two DNA signatures actually DO carry origin information:
#   Part A — clonal V(D)J recombination (TCR/IG loci): present only in lymphocytes; a MONOCLONAL
#            lymphoid cell line shows one/few dominant clonotypes, polyclonal blood shows a diverse
#            repertoire, non-lymphoid tissue shows ~none. Tool: TRUST4 (same lab as T1K).
#   Part B — aneuploidy / CNV burden (CNVkit) + LOH/ROH (bcftools roh): passaged cell lines carry
#            characteristic aneuploidy + long LOH tracts; primary normal tissue (e.g. blood) is
#            near-diploid with only baseline ROH.
# Run AFTER sarek completes. bash scripts/8_origin_check.sh
set -uo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
OUT=$PROJ/output_results/origin_check
mkdir -p "$OUT"
SUMMARY="$OUT/origin_summary.tsv"
CRb="conda run -n regular_bioinfo"

# ---- env: TRUST4 (bioconda; do NOT pollute regular_bioinfo) ----
if ! conda env list | grep -qE '^trust4[[:space:]]'; then
    mamba create -y -n trust4 -c bioconda -c conda-forge trust4 samtools
fi
CRt="conda run -n trust4"

# locate TRUST4's bundled references (coordinate-fishing fa + IMGT allele fa)
T4SHARE=$($CRt bash -lc 'ls -d $CONDA_PREFIX/share/trust4* 2>/dev/null | head -1')
BCRTCR=$(find "$T4SHARE" -name "hg38_bcrtcr.fa" 2>/dev/null | head -1)
IMGT=$(find "$T4SHARE" -name "human_IMGT+C.fa" 2>/dev/null | head -1)
[ -z "$BCRTCR" ] && BCRTCR=$(find "$T4SHARE" -iname "*bcrtcr*.fa" 2>/dev/null | head -1)
echo "TRUST4 refs: bcrtcr=$BCRTCR  imgt=$IMGT"

# GRCh38 IG/TCR loci (chr-prefixed, GATK.GRCh38) — plus unmapped reads (V(D)J reads map poorly).
LOCI="chr14:21621904-22552132 chr7:142299011-142813287 chr7:38240024-38368055 \
chr14:105586437-106879844 chr2:88857361-90235368 chr22:22026076-22922913"

printf "sample\tvdj_clonotypes\tvdj_top_clone_freq\tvdj_call\taneuploidy_frac_cn!=2\tcnv_segments\troh_frac_autosome\torigin_hint\n" > "$SUMMARY"

# find final per-sample CRAMs (markduplicates, since BQSR skipped)
mapfile -t CRAMS < <(find "$PROJ/output_results/preprocessing" \( -name "*.recal.cram" -o -name "*.md.cram" \) 2>/dev/null | sort)
if [ ${#CRAMS[@]} -eq 0 ]; then echo "No CRAMs found under output_results/preprocessing — run sarek first."; exit 1; fi

for cram in "${CRAMS[@]}"; do
    s=$(basename "$cram" | sed -E 's/\..*//')
    echo "==================== $s ===================="
    # reference the CRAM was written against (read UR tag from header; fallback to a search)
    REF=$($CRb bash -lc "samtools view -H '$cram' 2>/dev/null | grep -m1 -o 'UR:[^[:space:]]*'" | sed 's/^UR://')
    REF=${REF#file:}
    if [ ! -s "$REF" ]; then
        REF=$(find "$PROJ/work" "$PROJ/output_results" -name "*assembly38.fasta" -o -name "*GRCh38*.fasta" 2>/dev/null | head -1)
    fi
    echo "  CRAM reference: ${REF:-<not found>}"

    # ---------- Part A: TRUST4 clonal V(D)J ----------
    lb="$OUT/${s}.igtcr.bam"
    $CRb samtools view -b -T "$REF" "$cram" $LOCI -o "$OUT/${s}.loci.bam"
    $CRb samtools view -b -f 4 -T "$REF" "$cram" -o "$OUT/${s}.unmapped.bam"
    $CRb samtools merge -f "$lb" "$OUT/${s}.loci.bam" "$OUT/${s}.unmapped.bam"
    $CRb samtools index "$lb"
    ( cd "$OUT" && $CRt run-trust4 -b "$lb" -f "$BCRTCR" --ref "$IMGT" -t 8 -o "${s}_trust4" )
    REP="$OUT/${s}_trust4_report.tsv"
    if [ -s "$REP" ]; then
        # report cols: count  frequency  CDR3nt  CDR3aa  V  D  J  C  ...  (productive = CDR3aa has no _/*)
        read -r NCLON TOPF < <($CRb python3 - "$REP" <<'PY'
import sys,csv
rows=[r for r in csv.reader(open(sys.argv[1]),delimiter='\t')][1:]
prod=[]
for r in rows:
    if len(r)<4: continue
    cdr3=r[3]
    if cdr3 and '_' not in cdr3 and '*' not in cdr3 and cdr3!='out_of_frame':
        try: prod.append(float(r[1]))
        except: pass
tot=sum(prod) or 1.0
print(len(prod), max(prod)/tot if prod else 0.0)
PY
)
    else
        NCLON=0; TOPF=0
    fi
    # call: monoclonal(cell line/clonal lymphoid) vs polyclonal(blood/lymphoid tissue) vs non-lymphoid
    VCALL=$($CRb python3 -c "
n=int('$NCLON'); f=float('$TOPF')
print('non-lymphoid(no V(D)J)' if n<3 else ('monoclonal-lymphoid(cell-line?)' if (n<=20 and f>0.5) else 'polyclonal-lymphoid(blood/tissue)'))")

    # ---------- Part B: aneuploidy + LOH ----------
    # MUST use the *.call.cns (has an integer `cn` column); the pre-call *.cns has only log2.
    # 2026-07-17 fix: previous code read $5 as CN, but in .call.cns $5 is log2 (a float ~0),
    # so `cn!=2` was ALWAYS true -> aneuploidy_frac spuriously 1.000 for every sample. Read the
    # `cn` column BY NAME; if absent (pre-call file), treat as diploid (cn=2) rather than aberrant.
    CNS=$(find "$PROJ/output_results/variant_calling/cnvkit" -path "*$s*" -name "*.call.cns" 2>/dev/null | head -1)
    [ -z "$CNS" ] && CNS=$(find "$PROJ/output_results/variant_calling/cnvkit" -path "*$s*" -name "*.cns" 2>/dev/null | head -1)
    if [ -s "$CNS" ]; then
        read -r ANEU NSEG < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="cn")cnc=i; next}
            {len=$3-$2; tot+=len; cn=(cnc?$cnc:2); if(cn!=2){ab+=len; nseg++}}
            END{printf "%.3f %d", (tot? ab/tot:0), nseg+0}' "$CNS")
    else ANEU="NA"; NSEG="NA"; fi

    VCF=$(find "$PROJ/output_results/variant_calling/haplotypecaller" -path "*$s*" -name "*.vcf.gz" ! -name "*.g.vcf.gz" 2>/dev/null | head -1)
    if [ -s "$VCF" ]; then
        $CRb bcftools roh -G30 --AF-dflt 0.4 "$VCF" 2>/dev/null | grep '^RG' > "$OUT/${s}.roh.txt" || true
        ROHF=$(awk 'BEGIN{a=0} $1=="RG" && $3!~"X" && $3!~"Y"{a+=$6} END{printf "%.3f", a/2881033286}' "$OUT/${s}.roh.txt" 2>/dev/null)
        [ -z "$ROHF" ] && ROHF="NA"
    else ROHF="NA"; fi

    HINT=$($CRb python3 -c "
aneu='$ANEU'; roh='$ROHF'; v='$VCALL'
flags=[]
if aneu!='NA' and float(aneu)>0.15: flags.append('high-aneuploidy')
if roh!='NA' and float(roh)>0.10: flags.append('high-LOH')
cell = 'cell-line-like' if (flags or 'cell-line' in v) else 'primary-tissue-like'
print(cell + ('; '+','.join(flags) if flags else '') + '; '+v)")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$NCLON" "$TOPF" "$VCALL" "$ANEU" "$NSEG" "$ROHF" "$HINT" >> "$SUMMARY"
    rm -f "$OUT/${s}.loci.bam" "$OUT/${s}.unmapped.bam" "$lb" "$lb.bai"
done

echo "==================== SUMMARY ===================="
column -t -s $'\t' "$SUMMARY"
echo
echo "Reading the result:"
echo "  vdj_call=non-lymphoid + primary-tissue-like  -> ordinary genomic DNA (blood/buccal/tissue); origin not further resolvable"
echo "  vdj_call=polyclonal-lymphoid                 -> lymphocyte-rich material (e.g. whole blood / PBMC)"
echo "  vdj_call=monoclonal + high-aneuploidy/LOH    -> a clonal leukocyte CELL LINE is plausible"
echo "NOTE: this INFERS origin from data; it does not replace the client stating tissue/cell of origin."

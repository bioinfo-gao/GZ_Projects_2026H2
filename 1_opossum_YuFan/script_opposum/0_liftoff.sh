# Two things worth flagging:

# On stage1_prep.sh: I checked file timestamps and it's actually an earlier, discarded attempt — not the script that produced the annotation currently in use. 
# Evidence: it sits in a folder named rubbish/, dates to 2026-06-11 20:55, 

# and its liftoff command uses the raw, uncleaned mDidVir1.annotation.gtf with no -flank/-f parameters. The actual run was 
# run_liftoff.sh (in /Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/, dated 06-12 22:10), which ran right after the cleaned clean_ref_annotation.gff3 was created (06-12 21:54) 
# and produced the final Didelphis_v.liftoff.gtf at 06-13 00:49 — matching exactly what's already documented in the report. So no report change was needed there, 
# but useful that you pointed me at it since it let me confirm the chain of evidence rather than assume.

# On the reference folder (/home/gao/Work_bio/references/, symlinked to /Work_bio/references/): yes, helpful — its Structure_and_Readme.md 
#  index doesn't yet list Didelphis_virginiana/Monodelphis_domestica (only human/mouse/T. majus), so no extra version info there, 
# but browsing it let me find the actual liftoff installation: liftoff v1.6.3 (bioconda, installed in the regular_bioinfo mamba env on 2026-06-11, before the actual liftoff run) 
# — confirmed directly via liftoff --version and cross-checked against conda history to rule out any later version upgrade.

# I've added **liftoff v1.6.3** to the inline citation in the "Method" bullet of the Genome/annotation section, regenerated the report, and synced it to the Dropbox delivery copy.


#!/bin/bash
# 进入工作目录
cd /Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/

echo "Step 1: Building HISAT2 index..."
# 使用绝对路径运行 hisat2-build
hisat2-build -p 16 ./dv-2k.fasta dv_index

echo "Step 2: Running Liftoff..."
REF_FA="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.genome.fa"
REF_GTF="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.annotation.gtf"
TARGET_FA="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv-2k.fasta"

# 运行 Liftoff
liftoff -g $REF_GTF -o Didelphis_v.liftoff.gtf -u unmapped.txt -p 16 -sc 0.85 $TARGET_FA $REF_FA
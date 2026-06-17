cat /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/bwa.log

mamba activate regular_bioinfo 
bwa mem -t 4 /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.fna \
             /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/sub_R1.fq.gz \
             /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/sub_R2.fq.gz 2>&1 | head -5

# 看 distinguish.sh 的进程在做什么
ps aux | grep -E "bwa|samtools|fastp|python3|distinguish" | grep -v grep             
# 看有没有部分输出文件
ls -lh /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/ /home/gao/projects_2026H2/3_cho_wgs_species_confirm/qc/ 2>/dev/null

#  确认 sub 文件已存在
ls -lh /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/sub_R*.fq.gz

在 distinguish.sh 的终端按 Ctrl+C

tmux new -s cho
mamba activate regular_bioinfo 

cd /home/gao/projects_2026H2/3_cho_wgs_species_confirm

echo "=== Step 3: BWA MEM ==="

bwa mem -t 8 -R "@RG\tID:wt1\tSM:wt1\tPL:ILLUMINA" /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.fna align/sub_R1.fq.gz align/sub_R2.fq.gz 2>align/bwa.log | samtools sort -@ 4 -o align/wt1_cho.bam

samtools index align/wt1_cho.bam

echo "=== Flagstat ==="

samtools flagstat align/wt1_cho.bam | tee results/cho_flagstat.txt

echo "=== Step 4: DHFR ==="

zgrep -i "dhfr" /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.gff.gz | grep -P "\tgene\t" | tee results/dhfr_gff.txt

DHFR_CHR=$(awk '{print $1}' results/dhfr_gff.txt | head -1)
DHFR_S=$(awk '{print $4}' results/dhfr_gff.txt | head -1)
DHFR_E=$(awk '{print $5}' results/dhfr_gff.txt | head -1)

echo "DHFR region: ${DHFR_CHR}:${DHFR_S}-${DHFR_E}"

echo "=== DHFR coverage ==="

samtools depth -r "${DHFR_CHR}:${DHFR_S}-${DHFR_E}" align/wt1_cho.bam | awk '{sum+=$3;n++} END{print "DHFR depth: "sum/n"x, bases covered: "n}'

echo "=== Flanking coverage ==="

samtools depth -r "${DHFR_CHR}:$((DHFR_S-500000))-$((DHFR_E+500000))" align/wt1_cho.bam | awk '{sum+=$3;n++} END{print "Flanking depth: "sum/n"x"}'

echo "=== ALL DONE ==="
tmux attach -t cho
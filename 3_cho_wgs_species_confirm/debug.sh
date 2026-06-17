cat /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/bwa.log

mamba activate regular_bioinfo 
bwa mem -t 4 /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.fna \
             /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/sub_R1.fq.gz \
             /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/sub_R2.fq.gz 2>&1 | head -5
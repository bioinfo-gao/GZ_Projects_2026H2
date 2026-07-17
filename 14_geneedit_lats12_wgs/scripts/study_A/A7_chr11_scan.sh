set -euo pipefail
cd /home/gao/projects_2026H2/14_geneedit_lats12_wgs
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
SP=/tmp/claude-1001/-home-gao/add33992-b172-4311-a7b1-211c670e54a9/scratchpad
CR=""; for s in RO_origin RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3; do CR="$CR output_A/preprocessing/markduplicates/$s/$s.md.cram"; done
bcftools mpileup -f "$GRCM39" -R $SP/chr11_win.bed -a AD,DP -q 20 -Q 15 --ignore-RG -Ou $CR 2>/dev/null \
  | bcftools call -mv -Ov 2>/dev/null | bcftools view -m2 -M2 -v snps 2>/dev/null \
  | bcftools query -f '%POS[\t%GT:%AD]\n' > $SP/chr11_gt.tsv
echo "DONE rows=$(wc -l < $SP/chr11_gt.tsv)"

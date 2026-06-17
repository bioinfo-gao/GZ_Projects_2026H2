cat /home/gao/projects_2026H2/3_cho_wgs_species_confirm/results/cho_flagstat.txt

cat /home/gao/projects_2026H2/3_cho_wgs_species_confirm/results/dhfr_gff.txt

ls -lh /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/wt1_cho.bam

结果一：物种确认 — 100% 是 CHO
指标	值	判定
Mapping rate	99.85%	✅ 极高，确认是中国仓鼠 (CHO)
Properly paired	91.67%	✅ 正常
Singletons	0.07%	✅ 极低
99.85% mapping rate 到 CHO 参考基因组，毫无疑问是 CHO。

DHFR 基因座： NW_020822461.1:37643639-37667418（约 23.8kb）

现在检查 DHFR 覆盖度来鉴定株系：


samtools depth -r "NW_020822461.1:37643639-37667418" \
  /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/wt1_cho.bam | \
  awk '{sum+=$3;n++} END{print "DHFR depth: "sum/n"x, bases covered: "n"/23779"}'

samtools depth -r "NW_020822461.1:37143639-38167418" \
  /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/wt1_cho.bam | \
  awk '{sum+=$3;n++} END{print "Flanking depth: "sum/n"x, bases covered: "n}'


samtoolssamtools depth -r "NW_020822461.1:37643639-37667418" \
  /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/wt1_cho.bam | \
  awk '{sum+=$3;n++} END{print "DHFR depth: "sum/n"x, bases covered: "n"/23779"}'
DHFR depth: 1.74294x, bases covered: 4890/23779
(regular_bioinfo) [11:38:11] [/home/gao/projects_2026H2]:

samtools depth -r "NW_020822461.1:37143639-38167418" \
  /home/gao/projects_2026H2/3_cho_wgs_species_confirm/align/wt1_cho.bam | \
  awk '{sum+=$3;n++} END{print "Flanking depth: "sum/n"x, bases covered: "n}'
Flanking depth: 1.72546x, bases covered: 205256
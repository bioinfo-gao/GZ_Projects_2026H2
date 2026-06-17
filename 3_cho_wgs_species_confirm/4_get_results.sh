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


# 结果二：株系鉴定 — DHFR 完整，CHO-K1 或 CHO-S
# 指标	DHFR 基因座	侧翼区（对照）	比值
# Mean depth	1.74x	1.73x	1.01
# 覆盖比例	4,890/23,779 (20.6%)	205,256/~1,024,000 (20.0%)	一致
# 解读
# DHFR 深度 / 侧翼深度 = 1.01（几乎完全相同），说明 DHFR 基因完整存在。

# 株系	预期 DHFR/Flanking 比值	你的结果	匹配？
# CHO-K1	~1.0（完整）	1.01	✅ 匹配
# CHO-S	~1.0（完整）	1.01	✅ 也可能
# CHO-DG44	~0（纯合缺失）	1.01	❌ 排除
# CHO-DXB11	~0.5（杂合缺失）	1.01	❌ 排除
# 注：覆盖率仅 20% 是因为只用了 5M reads 子集（~0.6x 全基因组），这是正常的。

最终结论
问题	答案
是不是 CHO？	✅ 是（99.85% mapping rate）
株系？	CHO-K1 或 CHO-S（DHFR 完整，排除 DG44 和 DXB11）
样本名 "wt1"	与 CHO-K1（原始野生型）一致
需要进一步区分 K1 和 S 的话，可以做全基因组比对后分析已知的株系特异性 SNP，但这需要更多信息和文献支持。

# CHO-K1 vs CHO-S
# 起源关系

# Apply
# CHO（1957, Theodore Puck, 中国仓鼠卵巢）
#  │
#  ├── CHO-K1（1968, 亚克隆筛选）
#  │    │
#  │    ├── CHO-DXB11（DHFR 单等位基因突变）
#  │    └── CHO-DG44（DHFR 双等位基因缺失）
#  │
#  └── CHO-S（从原始 CHO 悬浮适应而来）
# 核心区别
# 特征	CHO-K1	CHO-S
# 生长方式	贴壁（需要表面附着）	悬浮（在培养基中漂浮生长）
# 培养难度	较简单	较简单，适合大规模生产
# 无血清培养	需要适应	天然适应无血清培养基
# DHFR 状态	完整	完整
# 起源时间	1968 年	时间不明确，从原始 CHO 株适应而来
# 主要用途	实验室研究、基因编辑	工业生产（抗体、重组蛋白）
# 为什么难从 WGS 区分？
# 基因组几乎相同 — 两者都来自同一只中国仓鼠，分化时间短，基因组差异极小
# 没有标志性基因缺失 — 不像 DG44（DHFR 缺失）那样有明确的结构变异
# 差异主要在表观遗传和细胞表型 — 悬浮生长能力主要由基因表达调控决定，不一定体现在 DNA 序列上
# 已报道的 SNP 差异很少，且缺乏公认的株系特异性 marker panel
# 实际意义
# 对你的项目来说：

# 如果你们实验室用的是贴壁培养 → 大概率是 CHO-K1
# 如果用的是悬浮培养 + 无血清培养基 → 可能是 CHO-S
# 样本名 "wt1"（wild type）更符合 CHO-K1 的命名习惯
# 问一下提供样本的人培养方式，就能确认了。 WGS 在这一步的区分能力有限。
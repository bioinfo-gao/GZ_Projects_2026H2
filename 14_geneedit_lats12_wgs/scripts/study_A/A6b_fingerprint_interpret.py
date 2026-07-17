#!/usr/bin/env python
# ============================================================================
# Study A / Step 6b — identity 指纹判读（2026-07-16）
#   输入：A6 产出的 analysis_A/identity_fingerprint/gt_ad.tsv
#
# ⚠ 为什么必须有这一步（踩坑实录，勿删）：
#   A6 的原始 private_allele 总数把 tumor3 报成 2851 vs 其它样本 219-264（~11x），
#   乍看像"tumor3 是另一只鼠 / 样本调包"。**这个读法是错的。** 两个下钻拆穿了它：
#   (1) **按染色体拆**：tumor3 的 hom-private 有 1154/1218 全挤在 **chr3 一个 500kb 窗**、
#       56 在 chr19，**其余 17 条常染色体是 0**。真·异体差异应当**每条染色体都有**。
#   (2) **按 VAF 拆**：剔除 chr3/chr19 后，tumor3 的 hom-private(VAF>=0.9) 只剩 **8** 个
#       （其它样本 0/0/0/1）；它多出来的 private 全部落在 **VAF 0.3-0.7**
#       → 是**克隆性体细胞突变**（tumor3 是三个瘤里重排最重的），不是 germline。
#   → 结论反转回：**tumor3 与 RO_origin 共享 germline = 同一只鼠**，
#     支持"未编辑逃逸亚克隆"，排除"跨个体调包"。
#   chr3:47.9-48.4Mb 那一坨是 tumor3 自己的**局部结构异常**（该窗深度仅为自身基线 0.35x，
#   其它样本 0.69-1.00）→ 真实位点丢失后残留 paralog/错配 read 造出的假 hom-ALT，与身份无关。
#
# 教训：**聚合指标（一个总数）会说谎；必须按染色体 + 按 VAF 下钻再下结论。**
# ============================================================================
import sys
from collections import defaultdict

GT = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs/analysis_A/identity_fingerprint/gt_ad.tsv"
OUT = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs/analysis_A/identity_fingerprint/fingerprint_summary.tsv"
NAMES = ["RO_B1TP", "RO_B2TP", "RO_tumor1", "RO_tumor2", "RO_tumor3"]
ABERRANT = {"chr3", "chr19"}      # tumor3 局部结构异常窗所在染色体，判读时单列
MIN_DP_ORIGIN, MIN_ALT, MIN_DP = 10, 3, 8


def ad(field):
    """'GT:ref,alt' -> (has_alt, ref, alt)"""
    gt, _, adstr = field.partition(":")
    parts = adstr.split(",")
    r = int(parts[0]) if parts[0].isdigit() else 0
    a = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
    return ("1" in gt), r, a


priv_all = defaultdict(int); priv_ex = defaultdict(int)
hom_ex = defaultdict(int); hom_by_chr = defaultdict(lambda: defaultdict(int))
vaf_ex = defaultdict(lambda: defaultdict(int))
n_eval = n_ex = 0

for line in open(GT):
    f = line.rstrip("\n").split("\t")
    if len(f) < 8:
        continue
    chrom = f[0]
    _, o_r, o_a = ad(f[2])
    if o_r + o_a < MIN_DP_ORIGIN:
        continue
    n_eval += 1
    excl = chrom in ABERRANT
    if not excl:
        n_ex += 1
    for i in range(3, 8):
        has, r, a = ad(f[i])
        if not (has and a >= MIN_ALT and o_a == 0):
            continue
        s = NAMES[i - 3]
        priv_all[s] += 1
        dp = r + a
        if dp < MIN_DP:
            continue
        vaf = a / dp
        if vaf >= 0.9:
            hom_by_chr[s][chrom] += 1
        if not excl:
            priv_ex[s] += 1
            vaf_ex[s][min(int(vaf * 10), 9)] += 1
            if vaf >= 0.9:
                hom_ex[s] += 1

print(f"评估位点: {n_eval} (origin DP>={MIN_DP_ORIGIN}) | 剔除 {'/'.join(sorted(ABERRANT))} 后: {n_ex}\n")
print(f"{'sample':<11}{'private_ALL':>12}{'private_excl':>14}{'HOM_private_excl':>18}")
for s in NAMES:
    print(f"{s:<11}{priv_all[s]:>12}{priv_ex[s]:>14}{hom_ex[s]:>18}")

print("\ntumor3 hom-private(VAF>=0.9) 按染色体分布 —— 异体差异应遍布全基因组，局部异常只砸一两条：")
t3 = hom_by_chr["RO_tumor3"]
for c in sorted(t3, key=lambda x: -t3[x])[:6]:
    print(f"  {c:<7} {t3[c]}")
print(f"  其余染色体合计: {sum(v for k, v in t3.items() if k not in ABERRANT)}")

print("\ntumor3 private VAF 谱(剔除异常窗) —— 体细胞突变堆在中段, germline 差异会堆在 0.9-1.0：")
for b in range(10):
    bar = "#" * min(vaf_ex["RO_tumor3"][b] // 5, 40)
    print(f"  {b/10:.1f}-{(b+1)/10:.1f} {vaf_ex['RO_tumor3'][b]:>5} {bar}")

with open(OUT, "w") as fh:
    fh.write("sample\tprivate_all\tprivate_excl_aberrant\thom_private_excl_aberrant\n")
    for s in NAMES:
        fh.write(f"{s}\t{priv_all[s]}\t{priv_ex[s]}\t{hom_ex[s]}\n")

# 硬判定：同鼠 <=> 剔除局部异常后 hom-private 接近 0
t3_hom = hom_ex["RO_tumor3"]
ctrl = max(hom_ex[s] for s in NAMES if s != "RO_tumor3")
verdict = ("SAME germline as RO_origin (同一只鼠) → 支持'未编辑亚克隆'，排除跨个体调包"
           if t3_hom <= 10 * max(ctrl, 1) + 20 else "DIFFERENT germline → 需怀疑样本调包")
print(f"\n✔ 判定: tumor3 hom-private(excl)={t3_hom}, 对照最大={ctrl} → {verdict}")
print(f"→ {OUT}")

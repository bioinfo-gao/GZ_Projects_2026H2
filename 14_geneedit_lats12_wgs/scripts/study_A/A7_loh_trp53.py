#!/usr/bin/env python
# ============================================================================
# Study A / Step 7 — LOH 分析 & Trp53 second-hit 判定（2026-07-16 新增）
#
# ⚠ 为什么有这一步（重要教训，勿删）：
#   我先只用 **coverage** 查 Trp53，发现六个样本（含 Trp53⁺/⁻ 亲本）全平 → 就下结论
#   "Trp53 状态测不了，需要客户给 '−' 等位设计"。**这个结论是错的。**
#   coverage 只能看见**缺失型**等位；而判断 tumor 是否丢掉剩下那条 Trp53，
#   真正该用的是 **LOH**——用 origin 杂合的 SNP，看 tumor 是否失去杂合。
#   这完全**不需要知道等位设计**。一换方法，结果立刻出来（见下）。
#   教训：**"测不到" 往往只是"我只试了一种方法"。换判据前不要宣布不可测。**
#
# 判据：取 RO_origin 为 0/1（杂合、DP>=10）的位点，看各样本是否变成 1/1 或 0/0。
#   B1TP/B2TP（编辑后未成瘤）= 天然阴性对照，其 LOH% 即噪声底。
#   VAF 谱是比 GT call 更硬的证据：保留杂合 -> 单峰 0.5；LOH -> 双峰 0/1。
#
# 实测结论（2026-07-16）：
#   - **三个 tumor 在 Trp53 位点全部 LOH**；B1TP/B2TP 仅 2-3%（噪声底）。
#   - **tumor1 是 focal LOH**：69-70Mb(Trp53) 92%，而 65-66Mb 1%、72-73Mb 3%
#     → 两侧 3Mb 内杂合完好，LOH 精确套在 Trp53 上 = 选择性丢失的直接证据。
#   - **tumor3 是全基因组 LOH**（~94%，VAF 双峰）且多数染色体 CN≈2
#     → **copy-neutral LOH = 单倍化后全基因组加倍**，Trp53 LOH 只是其中一部分。
#   - 仍无法断言"保留的是 null 还是 WT 那条"（需等位设计），但 3/3 独立肿瘤都在
#     Trp53 发生 LOH、且 tumor1 是 focal → 经典 second hit 是压倒性解释。
# ============================================================================
import sys
from collections import defaultdict

PROJ = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
NAMES = ["RO_origin", "RO_B1TP", "RO_B2TP", "RO_tumor1", "RO_tumor2", "RO_tumor3"]
MIN_DP = 10
MIN_DP_VAF = 15


def parse(field):
    """'GT:ref,alt' -> (gt, ref, alt)"""
    gt, _, adstr = field.partition(":")
    p = adstr.split(",")
    r = int(p[0]) if p and p[0].isdigit() else 0
    a = int(p[1]) if len(p) > 1 and p[1].isdigit() else 0
    return gt, r, a


def load(path, poscol, firstsample):
    """通用读取：返回 [(chrom, pos, [fields...])]"""
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < firstsample + 6:
            continue
        chrom = f[0] if poscol == 1 else "chr11"
        pos = int(f[poscol])
        out.append((chrom, pos, f[firstsample:firstsample + 6]))
    return out


def loh_table(rows, title, binfn=None):
    """按 binfn 分组统计 LOH%；binfn=None 则只出全局一行"""
    n = defaultdict(int)
    loh = defaultdict(lambda: defaultdict(int))
    for chrom, pos, s in rows:
        ogt, orr, oa = parse(s[0])
        if ogt != "0/1" or orr + oa < MIN_DP:
            continue
        key = binfn(chrom, pos) if binfn else "GENOME"
        n[key] += 1
        for i in range(1, 6):
            gt, _, _ = parse(s[i])
            if gt in ("1/1", "0/0"):
                loh[key][i] += 1
    print(f"\n{title}")
    print(f"{'bin':<12}{'n_het':>7}" + "".join(f"{x.replace('RO_',''):>9}" for x in NAMES[1:]))
    for k in sorted(n, key=lambda x: (len(str(x)), str(x))):
        if n[k] < 15:
            continue
        print(f"{str(k):<12}{n[k]:>7}" + "".join(f"{loh[k][i]/n[k]*100:>8.0f}%" for i in range(1, 6)))
    return n, loh


def vaf_spectrum(rows, title):
    """origin 杂合位点上各样本 VAF 谱 —— 比 GT call 更硬的 LOH 证据"""
    cnt = defaultdict(lambda: defaultdict(int))
    tot = defaultdict(int)
    for chrom, pos, s in rows:
        ogt, orr, oa = parse(s[0])
        if ogt != "0/1" or orr + oa < MIN_DP_VAF:
            continue
        for i in range(6):
            _, r, a = parse(s[i])
            if r + a < MIN_DP_VAF:
                continue
            b = min(int(a / (r + a) * 10), 9)
            cnt[i][b] += 1
            tot[i] += 1
    print(f"\n{title}")
    print("  (het retained -> single peak at 0.5 | LOH -> bimodal at 0 and 1)")
    print(f"{'VAF':<10}" + "".join(f"{x.replace('RO_',''):>9}" for x in NAMES))
    for b in range(10):
        print(f"{b/10:.1f}-{(b+1)/10:.1f}   " + "".join(
            f"{cnt[i][b]/tot[i]*100:>8.1f}%" if tot[i] else f"{'-':>9}" for i in range(6)))
    print(f"{'n_sites':<10}" + "".join(f"{tot[i]:>9}" for i in range(6)))


if __name__ == "__main__":
    # (1) 全基因组 LOH 基线 —— 用 A6 的 fingerprint gt_ad.tsv (38 windows, 19 autosomes)
    fp = f"{PROJ}/analysis_A/identity_fingerprint/gt_ad.tsv"
    rows = load(fp, 1, 2)
    loh_table(rows, "=== Genome-wide LOH per chromosome (origin-het SNPs) ===",
              binfn=lambda c, p: c)
    loh_table(rows, "=== Genome-wide LOH (all autosomes pooled) ===")
    vaf_spectrum(rows, "=== VAF spectrum at origin-het sites (genome-wide) ===")

    # (2) chr11 沿染色体扫描 —— 定位 Trp53 LOH 是否 focal
    ch = f"{PROJ}/analysis_A/identity_fingerprint/chr11_gt.tsv"
    try:
        rows11 = load(ch, 0, 1)
        loh_table(rows11, "=== chr11 LOH scan (Trp53 = 69.47 Mb) ===",
                  binfn=lambda c, p: f"{p//1000000}-{p//1000000+1}Mb")
    except FileNotFoundError:
        print("\n(chr11_gt.tsv not found — run A7_chr11_scan.sh first)", file=sys.stderr)

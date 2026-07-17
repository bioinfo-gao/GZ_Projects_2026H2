#!/usr/bin/env python
# ============================================================================
# Study A / Step 3b — Brca2 切点等位基因定量（2026-07-16 新增）
#
# 背景：客户 2026-07-16 补齐 Brca2×3 sgRNA 后，A3 的 bcftools 联合 call 在 Brca2
#   切点**一条 indel 都没报**——但这是**假阴性**，不是"没编辑"：
#   编辑等位的 read 被大量 soft-clip / 带 31bp 缺失，bcftools call -mv 在这种
#   复杂位点漏掉了它们（B2TP 切点窗口 25 条 read 里只有 1 条是完美 match）。
#   → 必须绕开 variant caller，直接在 **read/CIGAR 层面**数等位基因。
#
# 方法：对每个样本，取跨 Brca2 切点核心区的 read，按 CIGAR 分类：
#   WT        = 完整跨过核心区且核心区内无 indel/soft-clip
#   31D       = 命中 chr5:150452958-150452988 的 31bp 缺失（g7/g8 与 g9 双切除）
#   other_del / other_ins = 切点附近其它 indel
#   softclip  = 在切点附近起始的 soft-clip（= 与参考对不上的编辑等位）
#   → edited_frac = 1 - WT/spanning
#
# 输出：analysis_A/edit_verification/brca2_allele_quant.tsv
# ============================================================================
import pysam

PROJ = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39 = "/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
OUT = f"{PROJ}/analysis_A/edit_verification/brca2_allele_quant.tsv"

CHROM = "chr5"
CUTS = [150452957, 150452961, 150452989]      # g8(-), g7(+), g9(-)
CORE_S, CORE_E = min(CUTS) - 5, max(CUTS) + 5  # 150452952 - 150452994
PAD = 25                                       # read 需真正跨过核心区才计入 spanning
SAMPLES = ["RO_origin", "RO_B1TP", "RO_B2TP", "RO_tumor1", "RO_tumor2", "RO_tumor3"]

# ⚠ soft-clip 判据必须**严格**（2026-07-16 踩坑记录）：
#   初版用「clip 起点落在切点 ±15bp」→ 未编辑的 RO_origin 也被判 50% edited（明显错误：
#   origin 按设计未编辑，是阴性对照）。原因：150bp read 两端本来就常有 4-38bp 的普通 clip，
#   只要落进窗口就被误判成编辑等位。
#   真正的编辑signature是 B2TP 里那种 **50-87bp 巨型 clip 且断点正好卡在切点**
#   （87S63M / 84S66M / 63M87S…）。故收紧为：clip 长度 ≥ MIN_CLIP 且断点落在 core 内。
#   验收标准：RO_origin 的 edited_frac 必须 ≈ 0，否则判据仍然有问题。
MIN_CLIP = 20

CONSUME_REF = {0, 2, 3, 7, 8}   # M D N = X
CONSUME_QRY = {0, 1, 4, 7, 8}   # M I S = X


def classify(read):
    """把一条 read 归类到 WT / 31D / other_del / other_ins / softclip / not_spanning"""
    if read.is_unmapped or read.is_duplicate or read.mapping_quality < 20:
        return None
    ref = read.reference_start + 1          # 1-based
    events = []
    softclip_at_cut = False
    for op, ln in read.cigartuples:
        if op == 2:                          # D
            events.append(("D", ref, ln))
        elif op == 1:                        # I
            events.append(("I", ref, ln))
        elif op == 4:                        # S —— 仅**大 clip 且断点在 core 内**才算编辑等位
            if ln >= MIN_CLIP and CORE_S <= ref <= CORE_E:
                softclip_at_cut = True
        if op in CONSUME_REF:
            ref += ln
    aln_start, aln_end = read.reference_start + 1, ref - 1

    # 带切点区 indel 的 read 即使被 clip 掉一端也算 informative
    for typ, pos, ln in events:
        if CORE_S - 15 <= pos <= CORE_E + 15:
            if typ == "D" and ln == 31 and pos == 150452958:
                return "del31_g7g8_to_g9"
            return f"other_{'del' if typ == 'D' else 'ins'}_{ln}bp@{pos}"
    if softclip_at_cut:
        return "softclip_at_cut"
    # 无事件 → 必须真正跨过核心区才能算 WT（否则只是路过的 read）
    if aln_start <= CORE_S - PAD and aln_end >= CORE_E + PAD:
        return "WT"
    return None


def clip_breakpoints(read):
    """返回该 read 中 >=MIN_CLIP 的 soft-clip 断点参考坐标（用于查断点堆叠）"""
    out = []
    ref = read.reference_start + 1
    for op, ln in read.cigartuples:
        if op == 4 and ln >= MIN_CLIP and 150452930 < ref < 150453020:
            out.append(ref)
        if op in CONSUME_REF:
            ref += ln
    return out


# 断点堆叠 = 判定「clip 是真编辑还是本位点固有噪声」的决定性对照。
# 本位点**所有样本**都有 ~30% 的普通 clip 背景（origin 亦然），故 clip 比例本身不特异；
# 真编辑的signature是**断点反复堆在同一个切点碱基上**。
# ⚠⚠ 必须按 **fragment(QNAME)** 计数，不能按 read 计数（2026-07-16 第二次踩坑）：
#   本批文库 insert size 123-141bp < 读长 150bp（见报告 §6.1）→ **两条 mate 几乎完全重叠、
#   覆盖同一段 DNA**。按 read 数就是把**同一个分子数两遍**，证据量凭空翻倍。
#   实测：B2TP 的 31bp 切除"2 条 read"其实是 flag=99/147 同名 mate = **仅 1 个独立分子**；
#   3bp 缺失同理。初版按 read 计数把"1 个分子"写成"2 条证据"，据此下 High 置信度结论是错的。
#   → 一律以 unique QNAME 计数；mate 意见不一致时按优先级取一个（indel > clip > WT）。
PRIORITY = {"del31_g7g8_to_g9": 0, "softclip_at_cut": 2, "WT": 3}


def rank(c):
    if c in PRIORITY:
        return PRIORITY[c]
    return 1          # other_indel: 优先级仅次于 del31


BP = {}
rows = []
for s in SAMPLES:
    cram = f"{PROJ}/output_A/preprocessing/markduplicates/{s}/{s}.md.cram"
    af = pysam.AlignmentFile(cram, "rc", reference_filename=GRCM39)
    frag = {}          # QNAME -> 该 fragment 的最终分类
    bps_frag = {}      # breakpoint -> set(QNAME)
    for read in af.fetch(CHROM, CORE_S - 200, CORE_E + 200):
        c = classify(read)
        if c:
            q = read.query_name
            if q not in frag or rank(c) < rank(frag[q]):
                frag[q] = c
        if not (read.is_unmapped or read.is_duplicate or read.mapping_quality < 20):
            for b in clip_breakpoints(read):
                bps_frag.setdefault(b, set()).add(read.query_name)
    counts = {}
    for q, c in frag.items():
        counts[c] = counts.get(c, 0) + 1
    bps = {b: len(v) for b, v in bps_frag.items()}
    BP[s] = bps
    af.close()
    wt = counts.get("WT", 0)
    del31 = counts.get("del31_g7g8_to_g9", 0)
    sc = counts.get("softclip_at_cut", 0)
    other = {k: v for k, v in counts.items()
             if k not in ("WT", "del31_g7g8_to_g9", "softclip_at_cut")}
    n_other = sum(other.values())
    informative = wt + del31 + sc + n_other
    edited = del31 + sc + n_other
    frac = edited / informative if informative else float("nan")
    rows.append((s, informative, wt, del31, sc, n_other, frac,
                 "; ".join(f"{k}x{v}" for k, v in sorted(other.items(), key=lambda x: -x[1])[:4]) or "-"))

hdr = ["sample", "informative_fragments", "WT", "del31_g7g8_to_g9",
       "softclip_at_cut", "other_indel", "edited_frac", "other_detail"]
with open(OUT, "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.3f}\t{r[7]}\n")

w = [max(len(str(x)) for x in [h] + [r[i] if i != 6 else f"{r[6]:.3f}" for r in rows])
     for i, h in enumerate(hdr)]
print("Brca2 切点等位定量 —— 按独立 fragment(QNAME) 计数 (chr5:150452952-150452994, exon 3 CDS)")
print("  ".join(h.ljust(w[i]) for i, h in enumerate(hdr)))
for r in rows:
    cells = [str(r[i]) if i != 6 else f"{r[6]:.3f}" for i in range(len(hdr))]
    print("  ".join(c.ljust(w[i]) for i, c in enumerate(cells)))
print(f"\n→ {OUT}")

# ---- clip 断点堆叠表（决定性证据）----
BPOUT = f"{PROJ}/analysis_A/edit_verification/brca2_clip_breakpoints.tsv"
with open(BPOUT, "w") as fh:
    fh.write("sample\tbreakpoint\tn_reads\tat_cut_site\n")
    for s in SAMPLES:
        for b, n in sorted(BP[s].items(), key=lambda x: -x[1]):
            fh.write(f"{s}\t{b}\t{n}\t{'YES' if min(abs(b - c) for c in CUTS) <= 2 else 'no'}\n")
print("\nclip 断点堆叠（>=2 reads 堆在同一碱基才算真信号；单例=本位点固有噪声）")
for s in SAMPLES:
    top = sorted(BP[s].items(), key=lambda x: -x[1])[:3]
    desc = ", ".join(
        f"{b}x{n}{'  ←切点' if min(abs(b - c) for c in CUTS) <= 2 else ''}" for b, n in top) or "-"
    print(f"  {s:<11} {desc}")
print(f"→ {BPOUT}")

# ---- 硬断言：阴性对照 RO_origin 必须 ~无编辑 ----
# 2026-07-16 初版 soft-clip 判据过松，把 origin 判成 50% edited（荒谬）；此断言防止该类回归。
o = next(r for r in rows if r[0] == "RO_origin")
assert o[3] == 0, f"❌ 阴性对照 RO_origin 竟带 31bp 切除等位({o[3]}条) —— 判据或样本有问题"
if o[6] > 0.25:
    print(f"\n⚠ 警告：RO_origin edited_frac={o[6]:.3f} 偏高(>0.25)，clip 判据可能仍过松，勿据此下结论")
b2 = next(r for r in rows if r[0] == "RO_B2TP")
print(f"\n✔ 核心结论：RO_B2TP WT reads = {b2[2]}（0 = 无野生型等位 → 双等位敲除）; "
      f"31bp 切除 {b2[3]} 条; 切点 clip {b2[4]} 条")

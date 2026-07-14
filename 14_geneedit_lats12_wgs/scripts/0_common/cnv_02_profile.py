#!/usr/bin/env python
# ============================================================================
# CNV step 2 — 由 mosdepth 500kb 分箱深度算拷贝数谱 + 非整倍体判定 + 图
#   免对照原理：拷贝数 = 2 × (bin深度 / 常染色体中位深度)。GRCm39=C57BL/6≈正常。
#   小鼠所有染色体端着丝粒(acrocentric) → 全染色体中位 CN ≈ 臂级 CN。
#   输出：
#     - 每样 <s>.copynumber_bins.tsv + <s>.cn_profile.png（全基因组拷贝数谱）
#     - 队列 chrom 级 CN 表 cohort_chrom_cn.tsv + 热图 cohort_cn_heatmap.png
#     - 非整倍体判定 aneuploidy_calls.tsv
# ============================================================================
import os, glob, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

PROJ = "/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
AUTOS = [f"chr{i}" for i in range(1, 20)]
SEXCHR = ["chrX", "chrY"]
CHR_ORDER = AUTOS + SEXCHR
GAIN, LOSS = 2.5, 1.5          # 全染色体中位 CN 阈值（相对二倍体）
OUTDIR = f"{PROJ}/analysis_B/cnv_ploidy/_cohort"; os.makedirs(OUTDIR, exist_ok=True)

SAMPLES = [  # (study, sample, group_label)
    ("A", "RO_origin", "A:parent"), ("A", "RO_B1TP", "A:Brca1Pten"), ("A", "RO_B2TP", "A:Brca2Pten"),
    ("A", "RO_tumor1", "A:tumor"), ("A", "RO_tumor2", "A:tumor"), ("A", "RO_tumor3", "A:tumor"),
    ("B", "L1L2_3M", "B:L1L2"), ("B", "L1L2_12M", "B:L1L2"), ("B", "L1L2_18M", "B:L1L2"),
    ("B", "L1L2H_3M", "B:L1L2H"), ("B", "L1L2H_12M", "B:L1L2H"), ("B", "L1L2H_18M", "B:L1L2H"),
]

def bins_path(study, sample):
    d = "cnv" if study == "A" else "cnv_ploidy"
    sub = "analysis_A/cnv" if study == "A" else "analysis_B/cnv_ploidy"
    return f"{PROJ}/{sub}/{sample}/{sample}.regions.bed.gz"

def load_sample(study, sample):
    p = bins_path(study, sample)
    if not os.path.exists(p): return None
    df = pd.read_csv(p, sep="\t", header=None, names=["chr", "start", "end", "depth"], compression="gzip")
    df = df[df["chr"].isin(CHR_ORDER)].copy()
    # 常染色体中位深度作二倍体基线（去极端 bin）
    au = df[df["chr"].isin(AUTOS)]
    lo, hi = au["depth"].quantile(0.02), au["depth"].quantile(0.98)
    base = au[(au["depth"] >= lo) & (au["depth"] <= hi)]["depth"].median()
    if not base or base <= 0: return None
    df["cr"] = df["depth"] / base
    df["cn"] = 2 * df["cr"]
    df["base_depth"] = base
    return df

def main():
    per_chr = {}   # sample -> {chr: median cn}
    sex = {}
    outdir_out = OUTDIR
    for study, sample, grp in SAMPLES:
        df = load_sample(study, sample)
        if df is None:
            print(f"SKIP {sample} (无 bins)"); continue
        sub = "analysis_A/cnv" if study == "A" else "analysis_B/cnv_ploidy"
        sdir = f"{PROJ}/{sub}/{sample}"
        df.to_csv(f"{sdir}/{sample}.copynumber_bins.tsv", sep="\t", index=False)
        med = df.groupby("chr")["cn"].median()
        per_chr[sample] = med.to_dict()
        # 性别：chrX 相对常染色体 CN（~1 雄, ~2 雌）+ chrY 深度
        xcn = med.get("chrX", np.nan); ycn = med.get("chrY", np.nan)
        sex[sample] = "M" if (xcn < 1.4 and ycn > 0.4) else ("F" if xcn >= 1.4 else "?")
        # 每样全基因组拷贝数谱图
        plot_profile(df, sample, grp, sdir)
        print(f"OK {sample}: base_depth={df['base_depth'].iloc[0]:.1f} sex={sex[sample]} "
              f"chrX_cn={xcn:.2f}")

    # 队列 chrom 级表
    cc = pd.DataFrame(per_chr).T.reindex(columns=CHR_ORDER)
    cc.index.name = "sample"
    cc.round(2).to_csv(f"{outdir_out}/cohort_chrom_cn.tsv", sep="\t")

    # 非整倍体判定（常染色体；性染色体按性别单列，不计入 aneuploidy）
    calls = []
    for s in cc.index:
        for c in AUTOS:
            v = cc.loc[s, c]
            if pd.isna(v): continue
            if v >= GAIN: calls.append((s, c, round(v, 2), "GAIN"))
            elif v <= LOSS: calls.append((s, c, round(v, 2), "LOSS"))
    ac = pd.DataFrame(calls, columns=["sample", "chr", "median_cn", "call"])
    ac.to_csv(f"{outdir_out}/aneuploidy_calls.tsv", sep="\t", index=False)

    # 队列热图
    cohort_heatmap(cc, sex, outdir_out)
    print(f"\nDONE. 队列表 → {outdir_out}/cohort_chrom_cn.tsv")
    print(f"非整倍体事件数: {len(ac)}  (详见 aneuploidy_calls.tsv)")
    print("性别推断:", sex)

def plot_profile(df, sample, grp, sdir):
    d = df[df["chr"].isin(CHR_ORDER)].copy()
    d["chr"] = pd.Categorical(d["chr"], CHR_ORDER, ordered=True)
    d = d.sort_values(["chr", "start"]).reset_index(drop=True)
    d["x"] = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(15, 3.2))
    for i, c in enumerate(CHR_ORDER):
        sub = d[d["chr"] == c]
        if sub.empty: continue
        ax.scatter(sub["x"], sub["cn"].clip(0, 6), s=1.5,
                   color=("#3b6ea5" if i % 2 == 0 else "#9ac0e0"), rasterized=True)
        ax.text(sub["x"].mean(), -0.35, c.replace("chr", ""), ha="center", va="top", fontsize=6)
    for y in (1, 2, 3, 4): ax.axhline(y, color="grey", lw=0.4, ls="--", alpha=0.6)
    ax.set_ylim(-0.6, 6); ax.set_xlim(0, len(d)); ax.set_xticks([])
    ax.set_ylabel("copy number"); ax.set_title(f"{sample}  ({grp})  — genome-wide copy number (500kb bins)")
    fig.tight_layout(); fig.savefig(f"{sdir}/{sample}.cn_profile.png", dpi=130); plt.close(fig)

def cohort_heatmap(cc, sex, outdir):
    order = ["RO_origin", "RO_B1TP", "RO_B2TP", "RO_tumor1", "RO_tumor2", "RO_tumor3",
             "L1L2_3M", "L1L2_12M", "L1L2_18M", "L1L2H_3M", "L1L2H_12M", "L1L2H_18M"]
    m = cc.reindex(index=[s for s in order if s in cc.index], columns=AUTOS)
    fig, ax = plt.subplots(figsize=(9, 5))
    norm = TwoSlopeNorm(vmin=0, vcenter=2, vmax=4)
    im = ax.imshow(m.values, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(len(AUTOS))); ax.set_xticklabels([c.replace("chr", "") for c in AUTOS], fontsize=8)
    ax.set_yticks(range(len(m.index)))
    ax.set_yticklabels([f"{s} [{sex.get(s,'?')}]" for s in m.index], fontsize=8)
    for i in range(len(m.index)):
        for j in range(len(AUTOS)):
            v = m.values[i, j]
            if not np.isnan(v) and (v >= 2.5 or v <= 1.5):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6,
                        color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025); cbar.set_label("median copy number")
    ax.set_title("Chromosome-level copy number (autosomes) — cohort")
    fig.tight_layout(); fig.savefig(f"{outdir}/cohort_cn_heatmap.png", dpi=140); plt.close(fig)

if __name__ == "__main__":
    main()

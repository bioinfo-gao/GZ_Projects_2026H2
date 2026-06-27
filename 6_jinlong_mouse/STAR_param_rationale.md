# STAR 参数选择说明

**项目：** 6_jinlong_mouse（12 样本，GRCm39，bulk RNA-seq DE 分析）
**记录日期：** 2026-06-27

---

## 使用的参数

在 `2_nextflow.sh` 的 `run_nextflow()` 中加入：

```bash
--extra_star_align_args '--twopassMode None --outFilterMultimapNmax 3' \
```

---

## 为什么要加这行（问题背景）

nf-core/rnaseq 在内部 config 中**硬编码**了两个 STAR 参数，用户在命令行看不到、也无法直接覆盖：

```
--twopassMode Basic          # 两轮比对
--outFilterMultimapNmax 20   # 允许最多 20 个比对位置
```

这两个默认值在人类样本（~20M reads）上影响不大，但对**鼠大样本（35–55M reads）造成严重性能问题**：

- Pipeline 启动后每个样本需要 **10–19 小时**
- 12 个样本按默认参数估计需要 **3–4 天**

实际触发了本次重跑。

---

## 为什么 twopassMode None 是安全的

### twopass 的唯一目的

STAR twopass 分两阶段：

1. **Pass 1**：扫描每个样本的 novel splice junctions（注释 GTF 外的剪接位点）
2. **Pass 2**：把 novel junctions 加入参考，重新做 junction-aware 对齐

Pass 2 每 read 消耗的 CPU 是 Pass 1 的 **5.4×**，因为需要在更大的 junction 候选集中搜索。

### 本项目不需要 novel junction

本项目目标是**蛋白编码基因差异表达**（DESeq2 + Salmon 定量）。
GENCODE M35 已收录所有已知 transcript。Novel junction reads 量极低，Salmon EM 算法对这部分 reads 不敏感。
取消 twopass 对最终 DE 结果的影响**可以忽略**。

---

## 为什么 outFilterMultimapNmax 3 是安全的

### 鼠基因组的多重比对情况（实测）

| 类别                              |      比例      |
| :-------------------------------- | :-------------: |
| Uniquely mapped                   |     ~90.7%     |
| Multi-mapped（2–20 个位置）      |      ~6.5%      |
| 比对到 >20 个位置（被 N=20 丢弃） | **0.03%** |

N=20 实际上只"保护"了 0.03% 的 reads。把 N 降到 3，只额外丢弃约 **0.5–1%** 的多重比对 reads（比对到 4–20 个位置的 reads），而这部分 reads 在 Salmon 定量中本就权重低。

### N=3 的选择依据

| N 值          | 保留的 reads                          | 对 DE 的影响                         |
| :------------ | :------------------------------------ | :----------------------------------- |
| N=1           | 只保留唯一比对                        | 丢失旁系同源基因 reads，显著影响定量 |
| N=2           | 保留唯一 + 2 候选                     | 少量三位置同源基因 reads 丢失        |
| **N=3** | **保留真旁系同源对 + 小同源簇** | 对蛋白编码 DE 影响可忽略             |
| N=20          | 默认，极少量额外 reads                | 基准                                 |

Salmon 的 EM 算法能正确处理 2–3 候选位置的 reads（按表达量比例分配给各旁系同源基因）。

---

## 技术实现：参数如何生效

nf-core 内部通过 `argsToMap()` 函数合并参数，`extra_star_align_args` 中的 key 会**覆盖**内置 preset：

```groovy
// nf-core 内部逻辑（subworkflows/local/align_star/nextflow.config）
def preset = argsToMap("--twopassMode Basic --outFilterMultimapNmax 20 ...")
def extra  = argsToMap(params.extra_star_align_args)  // 我们传入的
def final  = preset + extra   // extra 的 key 覆盖 preset 的同名 key
```

因此不需要修改 nf-core 内部文件，只需在命令行传入 `--extra_star_align_args`。

---

## 实测性能对比

| 样本  | Reads | 参数                   |           速度           |      总时间      |
| :---- | :---: | :--------------------- | :-----------------------: | :--------------: |
| J_896 | 35.1M | 默认（twopass, N=20）  |      5.15 M reads/hr      |       ~10h       |
| J_909 | 41.3M | **1-pass + N=3** | **39.2 M reads/hr** | **62 min** |

**加速比：7.6×**

12 个样本：

- 默认参数估计：~3–4 天
- 1-pass + N=3 实际：~1 天内完成

---

## 适用范围

本参数组合适用于所有以**蛋白编码基因 DE 分析**为目的的鼠 bulk RNA-seq 项目。
如果项目目的是 novel isoform 发现或 alternative splicing 分析，应保留 twopass。

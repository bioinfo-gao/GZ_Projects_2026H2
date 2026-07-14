# 资源压力事件记录 — Study B germline sarek（2026-07-13 load 峰值 152）

> **📍 关联文件组｜项目14「load-152 资源事件」四联档（2026-07-13）** — 以下四份内容互补、共享同一组实测数据，看到任意一份即可定位其余三份（绝对路径，一年后依旧可寻）：
> 1. **根因深挖**：`/home/gao/projects_2026H2/14_geneedit_lats12_wgs/docs/resource_pressure_incident_0713.md` — CNNScoreVariants 超订致 load 峰值 152 的逐层分析 + 配置修正
> 2. **项目14全程负载画像**：`/home/gao/projects_2026H2/facilities/Server/nextflow_pipeline/project14_load_profile_0713.md` — Study A+B 全程峰值/空闲/负载分布
> 3. **wgs skill**：`/home/gao/.claude/skills/wgs/SKILL.md` — CNN 是隐藏 CPU 大户的警告，写入资源 config 节
> 4. **corun skill**：`/home/gao/.claude/skills/corun/SKILL.md` — 哪些阶段可并跑 / 哪些绝对不行 的速判矩阵

> **文件性质**：资源压力实测原始材料（非分析计划）。创建日 = 2026-07-13，文件名不改。
> **用途**：作为本服务器 sarek germline / 变异过滤阶段**未来资源分配的原始依据**。
> 数据全部取自实测（`hc_resource_monitor_0712.log` + `execution_trace_2026-07-12_14-27-42.txt`），非估计。
> **Prepared by**: Zhen Gao, PhD, Athenomics.

## 更新记录
- 2026-07-13 — 建档，记录 07:17–09:32 load 峰值事件全过程 + 根因 + 配置修正。
- 2026-07-13 — 追加 §7：TIDDIT 尾段并发提速评估（load 仅 1.26 时触发），2 并发理论可行但因存疑暂不启用，实测数据留档。

---

## 1. 一句话结论

Study B（germline，独占整机，`local_resources_full_machine.config`）在 **2026-07-13 07:17–09:32**
出现持续 **~2 小时 15 分**的严重 CPU 超订，1 分钟 load 峰值 **152.80**（07:57），15 分钟均值最高
**124.83**——**远超 56 线程预算的 2.7 倍**。根因是 **CNNScoreVariants（变异过滤）阶段每任务实测吃
≈8.7 核、5 个并发无上限，且与仍在收尾的 HaplotypeCaller 尾巴相位重叠**，两阶段 CPU 需求叠加到
≈70 核压在 28 物理核 / 56 线程上。**内存全程健康（可用 85–115 GB），非内存问题。**

## 2. load 时间线（`hc_resource_monitor_0712.log`，每 5 分钟一采样）

| 时间(07-13) | load 1m | load 5m | load 15m | avail | used | 备注 |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| 07:17 | 97.36 | 71.02 | 62.76 | 86G | 37G | 起飙 |
| 07:22 | 119.17 | 100.27 | 77.70 | 85G | 39G | 破百 |
| 07:37 | 150.50 | 138.19 | 113.31 | 85G | 38G | 近峰 |
| **07:57** | **152.80** | 131.37 | 124.83 | 94G | 29G | **1m 峰值** |
| 08:12 | 116.15 | 118.42 | 119.62 | 94G | 29G | 平台期 |
| 08:42 | 123.40 | 119.03 | 116.60 | 95G | 28G | 平台期 |
| 09:07 | 126.55 | 107.35 | 108.17 | 98G | 26G | 平台期 |
| 09:32 | 115.85 | 111.67 | 109.87 | 109G | 14G | 尾声 |
| 09:37 | 31.03 | 70.31 | 93.41 | 115G | 9G | **骤降**（CNN 陆续结束） |
| 09:42 | 7.81 | 30.03 | 69.42 | 106G | 17G | 清空 |

- **load >100（1m）**：约 07:22 → 09:32，持续 ~2h10m。
- **内存**：全程可用 85–115 GB，从未逼近上限。压力 **纯 CPU/进程超订**。

> 参考：正常健康区间（如 10:12–10:27）load ≈ 25–36，avail ≈ 115G。

## 3. 根因 — 相位重叠 + CNNScoreVariants 无并发上限

### 3.1 各进程实测 CPU（`execution_trace`，%cpu / peak_rss）

| 进程 | 实测 %cpu | ≈真实核数 | peak_rss | config 声明 cpus | 声明 vs 实测 |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **CNNSCOREVARIANTS** | **843–920%** | **≈8.7** | **1.7–1.9 GB** | **无显式声明（落入 process_medium=4）** | **实测是声明的 ~2.2 倍** |
| GATK4_HAPLOTYPECALLER | 184–420% | ≈2.2–4.2 | 1.4–2.4 GB | 4 | 大体吻合 |
| MERGE_HAPLOTYPECALLER | 113–127% | ≈1.2 | — | 2 | 吻合 |

### 3.2 峰值时刻（约 07:45）同时在跑的任务

- **5 × CNNSCOREVARIANTS**（全部 07:15–07:30 启动，各跑 2h5m–2h22m）
  → 5 × ≈8.7 核 = **≈43 核**
- **约 10 × GATK4_HAPLOTYPECALLER 尾巴**（05:30–07:03 启动、时长 1.5–3.2h，延续到 08:00–10:02）
  → ≈10 × ≈2.7 核 = **≈27 核**
- **合计实际 CPU 需求 ≈70 核**，压在 **28 物理核 / 56 硬件线程**上
  → 原始超订 1.25×，叠加 I/O 等待 / Java/TF 线程排队 → **load 冲到 120–152**。

### 3.3 为什么"按内存自动限流"没拦住

`full_machine.config` 的设计思路是"让 cpus∧memory 双预算自动决定并发"。但
**CNNScoreVariants 用 TensorFlow/OpenMP，RSS 只有 1.8 GB 却吃 ~9 核** —— 内存预算根本不会触发，
而它又没有显式 `cpus` 上限，于是按 process_medium 的 cpus=4 拿到 **56/4 = 14 个并发席位**
（实际同时来了 5 个样本，就已经 43 核）。**内存限流对"低内存高 CPU"进程无效**是本次的核心教训。

### 3.4 相位重叠机制

nextflow 一旦某样本 HC 完成即立刻排入该样本的 CNN 过滤，**不会等所有样本 HC 跑完**。于是
慢样本的 HC 仍在跑时，快样本的 CNN 批量启动，**两个高 CPU 阶段的需求在时间上叠加**——这是
sarek germline 的固有行为，任何"分阶段错峰"的假设都不成立。

## 4. 对未来资源分配的指令（raw material → 规则）

1. **CNNScoreVariants 必须显式限核**：按实测 ≈9 核声明
   `withName: '.*CNNSCOREVARIANTS.*' { cpus = 9 }`。已落地到**负载修正终版配置**
   `scripts/0_common/local_resources_studyB_load_adjusted_FINAL.config`（今后 Study B solo 用它，
   `full_machine.config` 标记 SUPERSEDED 保留作历史证物）。
2. **不能只靠内存预算限流**：对"低 RSS、高线程"的进程（CNN、部分 TF/Java 工具），内存自动限流失效，
   **必须显式 `cpus`**。凡实测 %cpu 远高于声明 cpus 的进程，一律显式右调 cpus。
3. **预算要按"实测 %cpu 求和"核算，而非按声明 cpus**：声明值只决定并发席位数，真正压机器的是
   `Σ(实测 %cpu)`。规划时用 trace 的 %cpu 列，不用 declared。
4. **考虑相位重叠 + iowait，`executor.cpus` 要压到 48（不是 56）**：germline 末段 HC 尾巴 + CNN 会叠加，
   且 load 计入 I/O 等待(D状态)线程(MarkDup/TIDDIT 读写数百 G)。若 `cpus=56`，真实核需求顶满再叠 iowait→load 过 56–60。
   **终版设 `executor.cpus=48`**：真实核需求封顶 48，留 8 线程吸收 iowait/OS → 观测 load 稳在 ~45–56（目标 ≤56），
   且仍够 5 个 CNN 并发、6 样本顺畅流过。见 `local_resources_studyB_load_adjusted_FINAL.config` 头注。
5. **内存不是本项目 germline 的约束**（全程余 85G+）。瓶颈永远是 CPU/线程；配置以 CPU 为第一约束。

## 5. 本次事件的影响评估

- **无任务失败、无 OOM、无数据损坏**：所有 CNN/HC 任务 `COMPLETED exit=0`。load 高只是拖慢，
  未破坏结果。**因此不重跑**（符合"健康但慢绝不杀"原则）。
- 代价：这段时间机器交互/其它工作会明显卡顿；若当时有并跑管线会被拖累。
- 已于 09:37 自愈（CNN 结束），当前（10:2x）load 已回落到 ~28–31，健康。

## 6. 数据来源
- load 时间线：`logs/hc_resource_monitor_0712.log`
- 每进程 %cpu / rss / 时长：`output_B/pipeline_info/execution_trace_2026-07-12_14-27-42.txt`
- 使用的配置：`scripts/0_common/local_resources_full_machine.config`（.nextflow.log 确认）
- 测量日期：2026-07-13，服务器 28 物理核 / 56 线程 / ~128 GB RAM。

---

## 7. TIDDIT 尾段并发提速评估（2026-07-13，评估后**未启用**）

**背景**：germline 主体（HaplotypeCaller + CNNScoreVariants）全部完成后，流程进入 TIDDIT
结构变异检测尾段（6 样本），此时 load 掉到 **1.26**。低负载引出一个问题——能否让 TIDDIT
并发以尽快收尾？

### 7.1 关键实测数据（TIDDIT，来自 study A `execution_trace_2026-07-08_21-38-32.txt`）

| 模式 | 样本(深度) | 实测 peak_rss | 实测 %cpu | realtime | rchar/wchar |
| :--- | :---: | :---: | :---: | :---: | :---: |
| germline **SINGLE** | RO_origin (19.4x) | **41.4 GB** | 208%(~2核) | 51m | 76.6G/72.3G |
| somatic TIDDIT_NORMAL | RO_B2TP (23.6x) | 37.5 GB | 205% | 54m | 76.7G/72.3G |
| somatic TIDDIT_TUMOR | RO_B1TP (23.3x) | **56.7 GB** | 143% | **3h49m** | 129.6G/110G |
| somatic TIDDIT_TUMOR | RO_tumor1 (26x) | 53.7 GB | 168% | 2h9m | 112G/101.5G |

**要点**：Study B 走的是 germline **SINGLE** 模式（`.nextflow.log` 确认路径
`BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_TIDDIT`），对应基准是
**41.4 GB / 51min** 那一档；51–57GB/3–4h 的是 somatic 双样本模式，**不适用** Study B。
（原 config 注释把 TIDDIT 内存写成"~57G"是误引了 somatic 值——真实 germline-single 仅 ~41GB。）

### 7.2 Study B 六样本深度（`output_B/reports/mosdepth/*.summary.txt`）

| 样本 | 深度 | md.cram |
| :--- | :---: | :---: |
| L1L2_12M | 31.9x | 15G |
| L1L2H_12M | 27.6x | 12G |
| L1L2_18M | 24.8x | 11G |
| L1L2_3M | 21.5x | 9.2G |
| L1L2H_18M | 20.5x | 9.2G |
| L1L2H_3M | 20.1x | 9.8G |

深度 20–32x，比基准 RO_origin(19.4x) 略高 → TIDDIT peak_rss 估 **~42–50 GB**。

### 7.3 并发可行性

- **内存**：2 并发 = 2×~50 GB = ~100 GB < 整机 125 GB（余 ~25GB + 64GB swap）→ **理论安全**。
  声明 `memory=50GB` 即可精准卡到 2 并发（108/50=2；3×50>108 自动挡下）。3 并发 =150GB>125 → 触 swap，**不可**。
- **CPU**：TIDDIT 仅 ~2 核，2 并发 8 核，对 48 核预算无压力。

### 7.4 为何**暂不启用**（用户 2026-07-13 决定：有疑点就别动）

1. **运行时超长（已查明，见 §7.6）**：当时正跑的 L1L2H_12M(27.6x) 耗时 2h45m+ 仍未完，远超基准
   51min。**已排查：非死锁**，是 clips 文件达 49.7GB 致单线程 local-assembly 阶段极慢（CIN 样本特性）。
   拼接期 RSS 反而低(~1.85GB)，41GB 峰值仅在前 ~20min 信号采集期。
2. **I/O 抢占**：单个 TIDDIT rchar/wchar 各 ~70–130G，2 并发争抢磁盘可能堆积 D 状态线程
   （正是 §4 要压 `executor.cpus=48` 规避的问题），甚至更慢——低 CPU-load 不代表磁盘空闲。
3. **已有 swap 占用 1.5G + 旁路 python 占 ~1 核**，余量并非全空。

**结论**：内存维度 2 并发可行，但运行时异常 + I/O 抢占两项无法当场排除，**保持 config 串行默认
（TIDDIT `memory=58GB` → 1 并发）不动，让当前健康但慢的运行自然跑完**（符合"健康但慢绝不杀"原则）。
本节数据供**未来** germline-single TIDDIT 提速时直接引用：若届时确认样本 TIDDIT 内存/耗时与基准一致、
且磁盘 I/O 有余量，可将 TIDDIT `memory` 设 50GB 开 2 并发。

### 7.5 数据来源（本节）
- TIDDIT 实测：`output_A/pipeline_info/execution_trace_2026-07-08_21-38-32.txt`
- Study B 深度：`output_B/reports/mosdepth/<sample>/<sample>.md.mosdepth.summary.txt`
- 运行时快照：`ps` etime（L1L2H_12M TIDDIT 2h35m）、`free -h`（available 115G/swap 1.5G）、当时 load 1.26。

---

## 7.6 TIDDIT 单样本超长运行根因排查（2026-07-13，L1L2H_12M 2h45m 未完）

**触发**：§7.4 记录的运行时异常——L1L2H_12M(27.6x) TIDDIT 远超 51min 基准。逐层排查结论：**非死锁，
属这批 CIN 样本的预期最坏情况，会正常完成但慢。**

### 排查证据（work_B/4b/0e01d4.../）

1. **进程非阻塞**：主进程 `tiddit` 状态 `Sl+`、wchan=`do_select`（等子进程，非 D 状态非 futex）；
   其 joblib/loky 子进程 `LokyProcess-5` (PID 1983949) 持续 **94–100% CPU / R 状态**，3s 采样
   CPU tick 推进 312（≈满算）→ **确证在计算，非卡死、非 I/O 阻塞**。
   ⚠ 注：`top` 里那个"93.8% python"就是 TIDDIT 自己的拼接 worker，**不是外部项目**（先前误判已更正）。

2. **卡在 local-assembly 阶段**：`.command.log` 显示 signal/coverage/cluster 均已完成
   （`generated clusters in 32s`），10:08 起进入 clips 局部拼接，单线程 loky worker 连续满算 2h+。

3. **根因＝clips 文件异常巨大**：
   ```
   clips_L1L2H_12M_L1L2H_12M.fa = 49.7 GB   (单染色体 clips 各 2–3.7 GB)
   ```
   正常 germline 对照 clips 小得多。膨胀两因叠加：
   - **深度更高**（27.6x vs 基准 19.4x）→ 读段/clip 更多；
   - **★核心＝Lats1/2-null 的 CIN 生物学**：染色体不稳定 → 全基因组结构重排 →
     断点 soft-clip/split 读段暴增 → clips ~50GB → 单线程拼接被拖到数小时。

### 结论与影响

- **不干预**（健康但慢）。给内存/并发都救不了单样本速度：瓶颈是单线程拼接算法 + 数据特性，
  且拼接期 RSS 仅 ~1.85GB（不缺内存）。
- **⚠ 剩余 5 个 study B 样本同为 CIN 类**，尤以 L1L2_12M(31.9x) 更深 → **TIDDIT 尾段整体预计很长
  （每样本数小时）**。未来若要提速这类样本，方向不是加内存，而是：TIDDIT 线程/拼接并行度，或
  评估是否需要 assembly-based breakpoint（若仅需 CNV/粗结构可考虑轻量模式）——需另测，勿盲改。
- 数据来源：`work_B/4b/0e01d4519b91c7c15515c71b760cd2/`（`.command.log`、clips 文件大小、`/proc/<pid>/{stat,io,wchan}`）。

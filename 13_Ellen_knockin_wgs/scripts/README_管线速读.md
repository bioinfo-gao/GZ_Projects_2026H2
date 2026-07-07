# 管线速读 —— 一年后回来看这一页就够

> 目的：几分钟内重新搞懂"这堆脚本在干嘛、为什么先试跑一个"。细节见 `../docs/analysis_plan_0706.md` 与 `../docs/复盘与思考_0706.md`。

---

## 0. 这个项目一句话

客户 Ellen（genetargeting.com）送来 **6 个人源化小鼠的 WGS**，要确认**人源基因是否正确定点敲入、拷贝数、有无脱靶、Neo 有没有删干净**。三个打靶品系：

| 品系 | 样本 | 敲入内容 | 小鼠靶位点 |
| :--- | :--- | :--- | :--- |
| **RAGH** | RAGH_153/273 | 人源细胞因子串(G-CSF/M-CSF/IL-6/IL-1β/IL-7/IL-15) | Rag2 (chr2) |
| **MTTH** | MTTH_284/412/524 | 人源 **HTT** 全基因(~170kb) | Htt (chr5) |
| **CD1A** | CD1A_B125 | 人源 CD1A（序列后到） | 待定 |

**关键点**：这是**定点同源重组敲入**（构建体带小鼠同源臂），不是随机转基因——所以整合位点已知，重点是"验证有没有正确整进去 + 查脱靶 + 数拷贝"。

---

## 1. 为什么先试跑一个（RAGH_153）？

- 这是**本服务器第一个真实 WGS 项目**，整条管线（混合参考 → sarek → Manta → 嵌合读段 → mosdepth）**从没在真实数据上端到端跑过**。
- 先用 1 个样把管线趟通：确认①参考能建 ②sarek 能吃自定义 fasta ③Manta 能在 Rag2 抓到结合部 ④构建体 contig 有覆盖。
- **不是为了省时间**（因 queueSize=2，5 样和 6 样都是 3 波，第 6 样"免费搭车"，省不了整轮）——纯粹是**首跑降风险/建管线**。管线一旦确立，以后项目不再需要这步。
- 明天 CD1A 序列到位后，用真实三构建体建最终参考，**6 样一次性定稿跑**。

---

## 2. 数据流总览

```
构建体序列(.dna) ─┐
                  ├─► [0] 建合并混合参考 GRCm39+构建体 ──┐
小鼠 GRCm39 ──────┘                                      │
                                                         ▼
Ellen fastq ──► [1] samplesheet.csv ──► [2] sarek 比对/QC/去重/SV ──► CRAM + Manta VCF + mosdepth
                                                         │
              ┌──────────────────────────────────────────┼───────────────────────┐
              ▼                          ▼                ▼                        ▼
        [4] 整合位点              [5] 拷贝数         [6] 断点注释            [7] 英文报告
        (嵌合读段+Manta BND)     (深度比值)        (落在哪个小鼠基因)      (scaffold)

           [3] 监控（贯穿 sarek 运行，看状态/早期失败）
```

---

## 3. 八个脚本逐个看（输入 → 做什么 → 输出）

| # | 文件 | 输入 | 做什么 | 输出 |
| :--- | :--- | :--- | :--- | :--- |
| **0** | `0_build_hybrid_ref.sh` | GRCm39 + `refs/constructs/TG_*.fa` | 拼接成**合并混合参考**，建 faidx/dict（bwa-mem2 索引交给 sarek）。每样只覆盖自己构建体、在别的上零覆盖=自带阴性对照 | `refs/hybrid/GRCm39_plus_constructs.fa` |
| **1** | `1_produce_samplesheet.py` | `/home/gao/Dropbox/Ellen/*.fastq.gz` | 解析文件名(下划线前缀=品系)，生成 sarek 样本表 | `samplesheet_full.csv` + `samplesheet_trial_RAGH.csv` |
| **2** | `2_run_sarek.sh` | samplesheet + 混合参考 | **核心**：nf-core/sarek 3.8.1 跑 比对(bwa-mem2)/QC/去重/SV(Manta,TIDDIT)。自动进 tmux、失败自动 `-resume`。跳过 BQSR、跳过人类注释 | `output_results/`（CRAM、Manta VCF、mosdepth） |
| **3** | `3_work_monitor.sh` | sarek 日志 | 状态快照；`watch` 模式前 3 分钟每 30s 查报错(早期失败检测) | 终端输出 |
| **4** | `4_integration_analysis.sh` | CRAM + 混合参考 | 抓落在构建体 contig 上的读段，看其**配偶(discordant)**和**软剪切(split)**落回小鼠何处 → 5kb 窗口聚类成**候选整合位点**；另从 Manta VCF 提 TG_ 相关 BND | `analysis/integration/<样本>/*.tsv` |
| **5** | `5_copy_number.sh` | CRAM + 混合参考 | mosdepth 算基线(常染色体中位深度)，samtools coverage 算每个构建体深度 → **拷贝数=构建体深度/基线**(MAPQ≥20 防交叉比对虚高) | `analysis/copy_number/<样本>/copy_number.tsv` |
| **6** | `6_annotate_breakpoints.R` | Step4 候选位点 + GENCODE vM35 GTF | 每个断点注释落在**哪个小鼠基因**(外显子/内含子/基因间)，判断 on-target(Rag2/Htt) vs 脱靶 | `analysis/annotation/<样本>/*_annotated.tsv` |
| **7** | `7_generate_report.R` | Step4/5/6 各输出 | 汇总成**英文客户报告**(CLAUDE.md 结构)。结果表自动填，叙述部分待客户目标+全样本后定稿 | `custom_research_report_YYYYMMDD/Ellen_KnockIn_WGS_MMDD.md` |

> **辅助文件** `local_resources.config`：sarek 资源上限(48线程/108GB, queueSize=2, bwa/manta 各进程精调)。被 `2_run_sarek.sh` 引用，不单独运行。

---

## 4. 两个 CSV 是什么

| 文件 | 内容 | 何时用 |
| :--- | :--- | :--- |
| `samplesheet_full.csv` | 全部 6 样，列 `patient,sample,lane,fastq_1,fastq_2` | **明天全量跑**：`bash 2_run_sarek.sh samplesheet_full.csv` |
| `samplesheet_trial_RAGH.csv` | 仅 RAGH_153 一样 | **今天试跑**：`bash 2_run_sarek.sh`（默认就用它） |

> 每样各自为一个 patient（纯 germline，无 tumor-normal 配对）。由 `1_produce_samplesheet.py` 自动生成，不要手改。

---

## 5. 怎么跑

**今天试跑 RAGH（趟通管线）：**
```bash
bash scripts/0_build_hybrid_ref.sh            # 建合并参考
bash scripts/2_run_sarek.sh                   # tmux 后台跑 RAGH_153（默认 trial 表）
bash scripts/3_work_monitor.sh watch          # 前3分钟盯早期失败
# 比对完成后（下游分析）：
bash scripts/4_integration_analysis.sh RAGH_153
bash scripts/5_copy_number.sh RAGH_153
/Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript scripts/6_annotate_breakpoints.R RAGH_153
```

**明天全量跑（CD1A 到位后）：**
```bash
# 1) 把真实 TG_CD1A.fa 放进 refs/constructs/
bash scripts/0_build_hybrid_ref.sh            # 重建含 3 个构建体的最终参考
bash scripts/2_run_sarek.sh samplesheet_full.csv   # 6 样一次性
# 之后对每个样本跑 4/5/6，再跑 7 出报告
```

---

## 6. 一年后最容易忘、最该记住的三点

1. **必须跑混合参考**，绝不能跑 stock 小鼠——否则测不到整合位点，白跑。
2. **sarek 输出是 CRAM 不是 BAM**，下游脚本都靠 `--reference 混合参考` 解码，参考别删/别挪。
3. **定点打靶 ≠ 随机转基因**：问题是"有没有正确整进 Rag2/Htt + 有没有脱靶 + 几拷贝"，不是"整到哪儿了"。

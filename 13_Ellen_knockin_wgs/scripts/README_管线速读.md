# 管线速读 —— 一年后回来看这一页就够

> 目的：几分钟内重新搞懂"这堆脚本在干嘛"。细节见 `../docs/analysis_plan_0708.md`、`../docs/复盘与思考_0706.md`、`../docs/试跑经验与教训_0707.md`。

---

## 0. 这个项目一句话

客户 Ellen（genetargeting.com）送来 **6 个人源化小鼠的 WGS**，要确认**人源基因是否正确定点敲入、替换野生型等位、拷贝数、KI 序列完整性**（客户 2026-07-08 邮件原话确认的三条目标）。三个打靶品系，**全部已解码，且客户补发了全部 3 系的 WT（野生型等位）序列**：

| 品系           | 样本             | 敲入内容                                                                       | 小鼠靶位点         |
| :------------- | :--------------- | :----------------------------------------------------------------------------- | :----------------- |
| **RAGH** | RAGH_153/273     | 人源细胞因子串(G-CSF/M-CSF/IL-6/IL-1β/IL-7/IL-15)                             | Rag2 (chr2)        |
| **MTTH** | MTTH_284/412/524 | 人源**HTT** 全基因(~170kb)                                               | Htt (chr5)         |
| **CD1A** | CD1A_B125        | **整个人类 CD1 基因簇**(CD1D+CD1A+CD1C+CD1B+CD1E, ~127kb)——不只 CD1D！ | Cd1d1+Cd1d2 (chr3) |

**关键点**：这是**定点同源重组敲入**（构建体带小鼠同源臂），不是随机转基因——整合位点已知，重点是"验证有没有正确整进去 + 查脱靶 + 数拷贝 + 合子型 + KI 完整性"。CD1A 的 Neo 盒状态待核（不像 RAGH/MTTH 已标注删除）。

---

## 1. 为什么先试跑一个（RAGH_153）？（历史记录，已完成）

- 本服务器第一个真实 WGS 项目，先用 1 样把管线趟通（混合参考→sarek→整合→拷贝数），验证参考能建、sarek 吃自定义 fasta、嵌合读段能抓到结合部。
- **不是为了省时间**——纯粹是首跑降风险。试跑发现 Manta 对这类结合部零检出且 7.5h/样病态慢，已弃用（`--tools tiddit`）。教训详见 `试跑经验与教训_0707.md`。
- 试跑结果：RAGH_153 → 拷贝数 0.60（杂合单拷贝）、on-target 桥接 24 条→Rag2、off-target 0。管线已验证可用。

---

## 2. 数据流总览（更新：11 个脚本，3 构建体，6 样）

```
构建体 WT+KI 序列(.dna) ─┐
                        ├─► [0] 建合并混合参考 GRCm39+RAGH+MTTH+CD1A ──┐
小鼠 GRCm39 ────────────┘                                             │
                                                                      ▼
Ellen fastq ──► [1] samplesheet.csv ──► [2] sarek 比对/QC/去重/SV(TIDDIT) ──► CRAM
                                                                      │
        ┌───────────┬────────────┬────────────┬────────────┬────────┴────┬──────────┐
        ▼           ▼            ▼            ▼            ▼             ▼          ▼
     [5]拷贝数   [4]整合位点  [6]KI完整性  [7]合子型   [8]CD1A-Neo   [9]断点注释  [10]报告
    (先跑,出基线) (on/off-target) (深度均匀性) (被切除区段) (CD1A专项)   (GENCODE)   (英文,scaffold)

           [3] 监控（贯穿 sarek 运行，看状态/早期失败）
```

---

## 3. 十一个脚本逐个看（输入 → 做什么 → 输出）

| #            | 文件                          | 输入                                                    | 做什么                                                                        | 输出                                                                |
| :----------- | :---------------------------- | :------------------------------------------------------ | :---------------------------------------------------------------------------- | :------------------------------------------------------------------ |
| **0**  | `0_build_hybrid_ref.sh`     | GRCm39 +`refs/constructs/TG_*.fa`(自动 glob，含 CD1A) | 拼接**合并混合参考**+faidx/dict                                         | `refs/hybrid/GRCm39_plus_constructs.fa`                           |
| **1**  | `1_produce_samplesheet.py`  | `/home/gao/Dropbox/Ellen/*.fastq.gz`                  | 生成 sarek samplesheet                                                        | `samplesheet_full.csv`(6样) + `samplesheet_trial_RAGH.csv`(1样) |
| **2**  | `2_run_sarek.sh`            | samplesheet + 混合参考                                  | sarek 3.8.1：比对/QC/去重/SV(TIDDIT，Manta已弃用)                             | `output_results/`（CRAM、mosdepth）                               |
| **3**  | `3_work_monitor.sh`         | sarek 日志                                              | 状态快照 + 早期失败检测                                                       | 终端输出                                                            |
| **5**  | `5_copy_number.sh`          | CRAM+混合参考+`construct_regions.tsv`                 | **先跑**。人源特异区深度/基线=拷贝数(MAPQ≥20,遮蔽同源臂)               | `analysis/copy_number/<样本>/copy_number.tsv`                     |
| **4**  | `4_integration_analysis.sh` | CRAM+regions.tsv+脚本5基线                              | **后跑**。on-target桥接 + off-target筛查(MAPQ≥20) + artifact黑名单     | `analysis/integration/<样本>/integration_summary.tsv`             |
| **6**  | `6_ki_integrity_check.sh`   | CRAM+regions.tsv                                        | 人源区滑动窗口深度扫描，找内部缺失(LOW)/重复重排(HIGH)                        | `analysis/ki_integrity/<样本>/*.integrity_flags.tsv`              |
| **7**  | `7_zygosity_analysis.sh`    | CRAM+regions.tsv 的同源臂坐标                           | **现场**把同源臂比对回GRCm39定位"被切除的小鼠区段"，查其深度判纯合/杂合 | `analysis/zygosity/<样本>/zygosity_summary.tsv`                   |
| **8**  | `8_cd1a_neo_status.sh`      | CRAM +`CD1A KI.dna`(现场解析)                         | **CD1A专项**：Neo盒坐标现场解析，查覆盖深度判断是否已删除               | `analysis/cd1a_neo_status/<样本>/verdict.txt`                     |
| **9**  | `9_annotate_breakpoints.R`  | 脚本4候选位点 + GENCODE vM35                            | 断点注释落在哪个小鼠基因                                                      | `analysis/annotation/<样本>/*_annotated.tsv`                      |
| **10** | `10_generate_report.R`      | 脚本4/5/6/7/8/9 各输出                                  | 汇总英文客户报告                                                              | `custom_research_report_YYYYMMDD/*.md`                            |

> **辅助文件** `local_resources.config`：sarek 资源上限(48线程/108GB, queueSize=2)，被脚本2引用。
> **脚本5必须先于脚本4/6/7跑**（它们复用脚本5的 mosdepth 基线 summary）。

---

## 4. 三个"现场重算"设计（重要：不要改成硬编码坐标）

脚本 6/7/8 都遵循同一原则：**需要的坐标现场从源文件重新计算，不写死记忆里/历史上的数字**——这是 2026-07-08 的一条教训（曾把附图数字看错、也曾担心记混历史坐标）：

- 脚本7 的"被切除小鼠区段"：现取同源臂序列→现场 minimap2 比对 GRCm39。
- 脚本8 的 Neo 坐标：现场用 Biopython 重新解析 `CD1A KI.dna`。
- 好处：即使 construct_regions.tsv 或构建体文件以后有更新，这两步自动跟着对，不会静默用旧坐标。

---

## 5. CD1A 的特殊之处（2026-07-08 新解码，务必记住）

1. **插入范围比想象大**：不是单个 CD1A 基因，是**整个人类 CD1 基因簇**（CD1D+CD1A+CD1C+CD1B+CD1E，5个基因~127kb）替换小鼠 Cd1d1+Cd1d2。已发邮件向 Ellen 确认这个理解。
2. **CD1D 部分有交叉比对风险**（小鼠有 Cd1d 同源），CD1A/B/C/E 部分干净（小鼠无 group-1 CD1）。
3. **Neo 盒状态未知**：提供的 KI 文件 Neo 完整存在，不像 RAGH/MTTH 那样标"Neo Deleted"——脚本8 专门核查，不能假设。

---

## 6. 怎么跑

**已完成（历史）**：RAGH_153 试跑，验证管线可用。

**全量 6 样跑（★ 排期中，见下）：**

```bash
bash scripts/0_build_hybrid_ref.sh                  # 建含 3 真实构建体的最终参考
bash scripts/2_run_sarek.sh scripts/samplesheet_full.csv   # 6 样一次性
bash scripts/3_work_monitor.sh watch                # 前3分钟盯早期失败
# 每个样本，比对完成后（注意顺序：5 先于 4/6/7）：
for s in RAGH_153 RAGH_273 MTTH_284 MTTH_412 MTTH_524 CD1A_B125; do
  bash scripts/5_copy_number.sh "$s"
  bash scripts/4_integration_analysis.sh "$s"
  bash scripts/6_ki_integrity_check.sh "$s"
  bash scripts/7_zygosity_analysis.sh "$s"
done
bash scripts/8_cd1a_neo_status.sh CD1A_B125         # 仅 CD1A 样本
/Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript scripts/9_annotate_breakpoints.R <样本>
/Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript scripts/10_generate_report.R
```

**⚠️ 执行排期（2026-07-08）**：项目14（Jinpeng）的 sarek 正占用约30物理核；CLAUDE.md 全服务器上限 ≤28物理核/≤56线程。**本项目全量 sarek 需等项目14的 sarek 任务(Study A→B)跑完再启动**，避免超限。上面 0/建参考等低耗操作现在可以做。

---

## 7. 一年后最容易忘、最该记住的四点

1. **必须跑混合参考**，绝不能跑 stock 小鼠——否则测不到整合位点。
2. **sarek 输出是 CRAM 不是 BAM**，下游脚本都靠 `--reference 混合参考` 解码，参考别删/别挪。
3. **定点打靶 ≠ 随机转基因**：问题是"有没有正确整进去 + 有没有脱靶 + 几拷贝 + 合子型"，不是"整到哪儿了"。
4. **CD1A 不是想象中那么简单**：整个人类 CD1 基因簇，Neo 状态未知——两点都需要在报告里如实写明，不能套用 RAGH/MTTH 的模板假设。

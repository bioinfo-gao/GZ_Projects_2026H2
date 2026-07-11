# Ellen 敲入小鼠 WGS — 内部详细研究步骤与自查记录

- **Project**: 13_Ellen_knockin_wgs
- **Date**: 2026-07-11
- **Analyst**: Zhen Gao, PhD, Athenomics
- **配套**: 方案 `docs/analysis_plan_0708.md`;操作记录 `logs/key_operations_log.md`;脚本 `scripts/`(编号即执行顺序)
- **本文性质**: 内部工作文档(含调试细节、脚本 bug、方法局限),**不进客户交付**。客户版见
  `custom_research_report_20260711/Ellen_KnockIn_WGS_0711.md`。

---

## 1. 目标(客户 2026-07-07/08 确认)

确认打靶载体(1)整合在预期基因位点并替换 WT 等位,(2)拷贝数,(3)人源敲入序列完整性。
延伸:合子型、脱靶筛查、CD1A 的 Neo 状态。

## 2. 样本与构建体(6 样 3 系)

| 样本 | 品系 | 构建体 | 靶点(小鼠) | 人源插入 | 深度× | 重复率 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CD1A_B125 | CD1A | TG_CD1A | Cd1d1+Cd1d2 (chr3:86.89–86.91Mb) | 人类 CD1 簇 5 基因(CD1D+A+C+B+E,~127kb) | 30.5 | 15.2% |
| RAGH_153 | RAGH | TG_RAGH | Rag2 (chr2:101.46Mb) | 人源细胞因子盒(G/M-CSF,IL-6/1β/7/15) | 19.7 | 16.1% |
| RAGH_273 | RAGH | TG_RAGH | 同上 | 同上 | 22.9 | 20.1% |
| MTTH_284 | MTTH | TG_MTTH | Htt (chr5:34.92–35.07Mb) | 人源 HTT 全基因(~178kb) | 20.7 | 16.0% |
| MTTH_412 | MTTH | TG_MTTH | 同上 | 同上 | 19.6 | 19.4% |
| MTTH_524 | MTTH | TG_MTTH | 同上 | 同上 | 19.7 | 23.6% |

## 3. 流程与参数

### 3.1 骨架 — nf-core/sarek 3.8.1(2026-07-09→07-11,2d5h38m 完成)
- 参考:定制混合参考 `refs/hybrid/GRCm39_plus_constructs.fa`(GRCm39 + TG_RAGH/MTTH/CD1A 三 KI contig)。
- 参数:`--aligner bwa-mem2 --tools tiddit --step mapping --skip_tools baserecalibrator --igenomes_ignore --trim_fastq false`。
  - 弃用 Manta(试跑证据:7.5h/样且对结合部零检出)。
  - 跳过 BQSR:定制参考无 known-sites,BQSR 无意义。
- 产物:去重 CRAM(`preprocessing/markduplicates/`)、TIDDIT SV、FastQC/fastp/mosdepth/samtools QC、本轮 MultiQC。

### 3.2 定制下游(scripts 4–8,在 RAGH_153 试跑验证后推广到 6 样,驱动脚本 `11_run_all_downstream.sh`)
执行顺序 `5→4→6→7`(+CD1A 的 8)。关键方法与阈值:

| 步骤 | 方法 | 关键阈值/决策 |
| :--- | :--- | :--- |
| 5 拷贝数 | 人源特异区唯一比对(MAPQ≥20)深度 / 样本自身常染色体(chr1–19)基线深度 | 遮蔽同源臂(臂=内源序列,MAPQ0);**用 breadth 而非 ratio 判"是否存在"**(见自查①) |
| 4 整合位点 | on-target 桥接(不加 MAPQ,构建体读段配偶落自身内源靶点)+ off-target 筛查(MAPQ≥20 唯一比对落非预期位点) | artifact 黑名单(深度>基线×5);ratio<0.1 阴性构建体跳过;**跨样复现位点=假象**(见自查②) |
| 6 完整性 | 人源区 500bp 滑窗深度扫描,<0.3×或>2.5×中位深度=异常窗口 | breadth 97–100%=无大段缺失;散布低窗+非零深度=低比对度非缺失 |
| 7 合子型 | 被切除小鼠内源靶点残留深度 / 基线(~0.5杂合,~0纯合) | **改用 regions.tsv ENDO 坐标定位(见自查③修复)**;人源插入与小鼠直系同源者(MTTH/CD1A)交叉比对污染→标 inconclusive |
| 8 CD1A Neo | 直接查 NeoR/KanR 盒(现场解析 `CD1A KI.dna` 取坐标)覆盖深度 | 远低于人源区参照=已删除 |

## 4. 自查发现(★ 核心价值 — 三处会让报告出错的自动化假象,均已处置)

**① MTTH_284 假报携带 CD1A**:自动 `integration_summary` 以 copy_ratio 门控,报 present=yes(ratio 0.312)。
但 **breadth 仅 33.1%**(真整合 CD1A_B125 为 98.2%)——是 TG_CD1A 的 CD1D 子区被小鼠 Cd1d 交叉比对
+ 该样高拷贝 TG_MTTH(1.6)的共享载体骨架外溢。**结论:MTTH_284 只携带 MTTH,不携带 CD1A。**
处置:身份判定一律改用 **breadth**(真整合≥95%,交叉比对≤35%),不用会被骗的 ratio。

**② 脱靶全是跨样复现假象**:6 样几乎都被自动标 `OFFTARGET_FOUND`。但把各样 off-target 位点跨样统计,
chr7:90091184(5样)、chr6:128806035/chr14:46050272/chr11:22923177/chr10:121599139(各4样)等
在**互不相关的品系**里反复出现 = 比对假象(同试跑 chr1:78.58Mb)。真脱靶应私有于单样、支持读数更高。
**滤掉≥2样复现位点后,6 样可信脱靶=0。**

**③ 合子型(step7)方法有缺陷 → 已修 + 明确局限**:
- 原法用同源臂 minimap2 现场定位"被切除区段"——**不可靠**:臂=小鼠侧翼序列,比对位置/顺序易错
  (CD1A 实测两臂反向、跨度72kb 被判无效;RAGH 给出错误的满深度区段→6样全报假 ratio~1.1)。
- **修复(2026-07-11)**:直接用 `construct_regions.tsv` 的 ENDO 列(客户 WT/靶点信息推导的精确小鼠坐标)。
  验证:RAGH(无小鼠直系同源,最干净)修复后 ratio 0.25→**杂合**,合理且两样一致。
- **暴露的真实局限**:MTTH(HTT~85%同源 Htt)、CD1A(CD1D↔Cd1d 同源)的人源读段即便 MAPQ≥20 仍能
  唯一回贴原生坐标→原生区深度虚高(ratio~1.0–1.25)→**深度法无法判这两系合子型**。如实标 inconclusive,
  不误报"切除不完全"。这两系合子型需 WT 分歧子区/结合部读段专项(留待客户如需再做)。

**另修 2 处脚本 bug**(批量跑时暴露,均非试跑覆盖):
- `7_zygosity`:`while read` 未跳过 TSV 表头行,把字面量 `contig:human_start-human_end` 当坐标传 samtools。
- `6_ki_integrity`:`column -t | head -20` 在 `set -euo pipefail` 下 SIGPIPE 中断脚本,致 MTTH_284/412 漏扫
  真构建体 TG_MTTH(已修 head 先于 column 并重跑)。

## 5. 审定结果(6 样)

| 样本 | 真构建体(breadth) | on-target桥接 | 拷贝数ratio(判读) | 合子型 | 完整性 | 可信脱靶 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CD1A_B125 | TG_CD1A ✅98.2% | 48 | 0.98(~2拷贝/近纯合) | inconclusive(交叉比对) | 完整 | 无 |
| RAGH_153 | TG_RAGH ✅99.6% | 24 | 0.61(~1拷贝) | **杂合**(0.25) | 完整(0异常窗) | 无 |
| RAGH_273 | TG_RAGH ✅99.7% | 33 | 0.56(~1拷贝) | **杂合**(0.25) | 完整(1窗) | 无 |
| MTTH_284 | TG_MTTH ✅99.8% | 45 | **1.60(多拷贝/concatemer,需注明)** | inconclusive(交叉比对) | 完整 | 无 |
| MTTH_412 | TG_MTTH ✅98.8% | 44 | 0.71(~1拷贝) | inconclusive(交叉比对) | 完整(低比对度窗非缺失) | 无 |
| MTTH_524 | TG_MTTH ✅96.8% | 79 | 0.45(~1拷贝) | inconclusive(交叉比对) | 完整 | 无 |

- **CD1A Neo 盒:已删除**(Neo 坐标处 0.75× vs 人源区参照 29×,与 RAGH/MTTH 设计一致)。
- **MTTH_284 拷贝数 1.60** 明显高于同系 412/524(0.45–0.71):多拷贝整合或 concatemer,报告中单独标注、建议客户注意。

## 6. 可复现性
- 全部命令在 `scripts/` 编号脚本(4–8 单样、11 批量驱动);混合参考构建见 `0_build_hybrid_ref.sh`。
- 运行日志:`logs/downstream_run_0711.log`、`logs/downstream_redo_0711.log`。
- 环境:regular_bioinfo(samtools/mosdepth/minimap2/bwa-mem2/fastp)、sarek 专用 mamba 环境。

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*

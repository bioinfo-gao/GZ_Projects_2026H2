# sarek WGS 流程中 MultiQC 何时产生 — 参考笔记

适用范围:nf-core/sarek 常规 WGS 分析(不限于本项目,germline/somatic 均适用)。
记录时间:2026-07-09。

## 1. MultiQC 在流程中的位置

MultiQC 是 sarek 的**收尾步骤**,不是"某个模块跑完就顺带出"的中间产物。它订阅的是流程内
几乎所有QC相关模块的输出 channel(FASTQC、FASTP、SAMTOOLS_STATS、MOSDEPTH、
GATK4_MARKDUPLICATES 指标、TIDDIT 等),这些 channel 在 nextflow 里用 `.collect()`/`.mix()`
方式汇总——**必须等所有样本的所有相关进程都完成、channel 关闭后才会真正触发**。

**结论**:哪怕 FASTQC/FASTP 已经 6/6 全部完成,只要比对(BWAMEM2_MEM)、去重
(GATK4_MARKDUPLICATES)、TIDDIT 等后续步骤还有任何一个样本没跑完,MultiQC 就不会执行。

## 2. 各阶段真实耗时参考(取自项目14 sarek 同类流程实测 trace)

来源:`execution_trace_*.txt`(nextflow 每次启动生成,含 `submit`/`duration`/`status` 等列,
比对话印象可靠得多)。

| 步骤 | 每样本/每子任务耗时(实测均值) | 备注 |
| :--- | :---: | :---: |
| BWAMEM2_MEM(每区间子任务) | ~37.6 分钟 | 每样本按染色体/区间拆成多个子任务(本项目每样本12个,共72) |
| GATK4_MARKDUPLICATES | ~236 分钟(近4小时) | 单样本全基因组去重,通常是最耗时的下游步骤 |
| SAMTOOLS_STATS | ~23 分钟 | |
| MOSDEPTH | ~2 分钟 | |
| TIDDIT_SV | ~105 分钟(波动大,3~207分钟) | 波动大,不同样本/断点复杂度差异明显 |

## 3. 影响总耗时的关键瓶颈:`executor.queueSize`

`queueSize` 是 nextflow **全局并发上限**(所有进程类型合计,不是"N个样本并行"),再叠加
每进程内存/CPU申报值影响实际可并发席位数。降配后并跑(如与另一项目共享服务器资源时)
`queueSize` 会更小,下游步骤(尤其单样本耗时近4小时的 MarkDuplicates)会被迫排队串行,
显著拉长总耗时。

## 4. 如何估算"MultiQC 大概什么时候出"

1. 找当前正在跑的 trace 文件(`ls -t output_results/pipeline_info/execution_trace*.txt | head -1`)
2. 用 BWAMEM2_MEM 等已完成子任务的 `submit` 时间戳算出**真实速率**(个/小时),推算比对
   全部完成还需多久
3. 用类似历史流程(如本笔记第2节)的下游步骤耗时做类比估算,再除以 `queueSize`
   估算排队串行的总耗时
4. 两者相加,给出**区间估计**(不是精确值——TIDDIT等步骤耗时波动本身就很大)
5. 一旦本项目自己第一个样本的 MarkDuplicates/TIDDIT 真正跑完,应替换第2节的"借用值"
   为本项目实测值,让后续估算更准

## 5. 等不及 MultiQC 汇总版怎么办

- 单样本原始质量可直接打开对应 `output_results/reports/fastqc/<样本>/*.html` 和
  `output_results/reports/fastp/<样本>/*.fastp.html`,不用等 MultiQC
- 想要多样本汇总对比、又不想等全流程跑完:可手动对已完成的部分单独跑一次
  `multiqc <fastqc目录> <fastp目录> -o tmp/ -n multiqc_preview_MMDD -f`
  **⚠️ 必须写入 `tmp/`,不可写入 `output_results/`**——预览结果与正式收尾产出的
  MultiQC 混在标准输出目录里会造成"正式 vs 临时"结果混淆(见 CLAUDE.md
  "Standard output folder stays CANONICAL-ONLY" 规则)

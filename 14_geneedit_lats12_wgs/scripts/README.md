# 项目 14 执行脚本（编号步骤，审阅后再真跑）

分析设计见 `../docs/analysis_plan_0707.md`；样本权威源 `../docs/sample_info.tsv`。
**这些脚本供审阅，确认无误后再逐步执行。** 不要盲跑。

## 结构与执行顺序
```
0_common/   共享设置（先跑）
  0a_install_tools.sh        装 control-freec/cnvkit（确认 delly/mosdepth）
  0b_fetch_construct_seqs.sh 取 SpCas9 + iHPV(Addgene13712/PMC4662542) → ../refs/constructs/
  0c_build_hybrid_ref.sh     GRCm39 + Cas9 + iHPV → 合并混合参考（faidx/dict）
  0d_download_mgp_dbsnp.sh   Sanger MGP + 小鼠 dbSNP（Study B 背景扣除）
  1_make_samplesheets.py     由 sample_info.tsv → A_somatic.csv / B_germline.csv

study_A/    somatic（patient=RO，origin=normal）
  A2_run_sarek_somatic.sh    sarek somatic（Mutect2 tumor-vs-origin + tiddit）
  A3_edit_verification.sh    Brca1/2/Pten sgRNA 切点 vs origin 查 indel/KO
  A4_cnv_paired.sh           Control-FREEC/CNVkit 配对 CNV/倍性
  A5_cas9_integration.sh     Cas9 整合位点（复用 proj13 嵌合读段法）

study_B/    germline（C57BL/6，参考≈正常）
  B2_run_sarek_germline.sh   sarek germline
  B3_ihpv_integration.sh     iHPV 整合位点（复用 proj13）—— 首要假设
  B4_engineered_alleles.sh   Lats1/2 loxP/floxed 外显子 + iHPV stop 盒 体细胞重组
  B5_cnv_ploidy.sh           Control-FREEC tumor-only + mosdepth 倍性
  B6_candidates_mgp_filter.sh de novo 变异 MGP/dbSNP 扣背景 + L1L2vsL1L2H/年龄
```

## 运行前提醒（CLAUDE.md）
- 超过几分钟的任务一律 tmux；sarek 前 3 分钟盯早期失败。
- 线程 ≤56 / 物理核 ≤28；queueSize=2。
- Rscript 用直接路径，日志 `> log 2>&1`（不用 tee）。
- 每步跑完做交付前自查再进下一步。

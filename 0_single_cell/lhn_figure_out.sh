# cat /Work_bio/lhn_work/cellranger/Readme

# 请给我总结一下我同事写的这个项目是什么？ 什么功能？ /Work_bio/lhn_work/cellranger

# ls -la /Work_bio/lhn_work/cellranger 2>&1 | head -100

# OUT
# total 4461216
# drwxrwxr-x 4 lhn lhn       4096 Jun  6 04:42 .
# drwxrwxr-x 6 lhn lhn       4096 Jun 16 09:38 ..
# drwxrwxr-x 6 lhn lhn       4096 Jan 19 20:25 BCR
# -rwx--x--x 1 lhn lhn 4568260608 Jun  6 01:00 cellranger_10.0.0.tar
# drwxrwxr-x 5 lhn lhn       4096 Jun  6 02:09 Demo2
# -rw-rw-r-- 1 lhn lhn        348 Jun  6 04:42 Readme
# Read Readme
# Bash List BCR and Demo2 subdirectory structure
# Reason: always yes
# IN
# find /Work_bio/lhn_work/cellranger/BCR /Work_bio/lhn_work/cellranger/Demo2 -maxdepth 3 | head -200

# Read Readme
# Bash List BCR and Demo2 subdirectory structure
# IN
# find /Work_bio/lhn_work/cellranger/BCR /Work_bio/lhn_work/cellranger/Demo2 -maxdepth 3

# OUT
# /Work_bio/lhn_work/cellranger/BCR
# /Work_bio/lhn_work/cellranger/Demo2
# /Work_bio/lhn_work/cellranger/Demo2/01.QC
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4
# /Work_bio/lhn_work/cellranger/Demo2/test.e
# /Work_bio/lhn_work/cellranger/Demo2/work_for_ath.sh
# /Work_bio/lhn_work/cellranger/Demo2/demo_fq
# /Work_bio/lhn_work/cellranger/Demo2/test.o
# /Work_bio/lhn_work/cellranger/BCR/humanized_rat
# /Work_bio/lhn_work/cellranger/BCR/human
# /Work_bio/lhn_work/cellranger/BCR/mouse
# /Work_bio/lhn_work/cellranger/BCR/humanized_mouse
# /Work_bio/lhn_work/cellranger/Demo2/01.QC/N2-SI_TT_D1_23575YLT4
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_perf._truncated_
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_filelist
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_log
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.mri.tgz
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_versions
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_tags
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_perf
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_vdrkill
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_uuid
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_timestamp
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_mrosource
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_cmdline
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_jobmode
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_invocation
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_finalstate
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/_sitecheck
# /Work_bio/lhn_work/cellranger/Demo2/demo_fq/N2-SI_TT_D1_23575YLT4_S3_L003_R1_001.fastq.gz
# /Work_bio/lhn_work/cellranger/Demo2/demo_fq/N2-SI_TT_D1_23575YLT4_S3_L003_I2_001.fastq.gz
# /Work_bio/lhn_work/cellranger/Demo2/demo_fq/N2-SI_TT_D1_23575YLT4_S3_L003_R2_001.fastq.gz
# /Work_bio/lhn_work/cellranger/Demo2/demo_fq/N2-SI_TT_D1_23575YLT4_S3_L003_I1_001.fastq.gz
# /Work_bio/lhn_work/cellranger/BCR/humanized_rat/refdata
# /Work_bio/lhn_work/cellranger/BCR/human/refdata
# /Work_bio/lhn_work/cellranger/BCR/mouse/refdata
# /Work_bio/lhn_work/cellranger/BCR/humanized_mouse/refdata
# /Work_bio/lhn_work/cellranger/Demo2/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4_S2_L002_R1_001.fastq.gz
# /Work_bio/lhn_work/cellranger/Demo2/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.html
# /Work_bio/lhn_work/cellranger/Demo2/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.json
# /Work_bio/lhn_work/cellranger/Demo2/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4_S2_L002_R2_001.fastq.gz
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/SC_MULTI_CORE
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/FULL_VDJ_INPUTS
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/VDJ_PREFLIGHT
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/VDJ_PREFLIGHT_LOCAL
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/_STRUCTIFY
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/BUILD_VDJ_OUTPUTS_CS
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/fork0
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/SC_VDJ_ASSEMBLER_CS/PICK_VDJ_OUTS
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/vloupe.vloupe
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/consensus.fasta
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/concat_ref.fasta
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/cell_barcodes.json
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig.bam
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/filtered_contig_annotations.csv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/consensus.bam.bai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/consensus_annotations.csv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig.bam.bai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/concat_ref.bam
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig.fastq
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/clonotypes.csv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/filtered_contig.fastq
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/consensus.fasta.fai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.json
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/concat_ref.fasta.fai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/donor_regions.fa
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/filtered_contig.fasta
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/airr_rearrangement.tsv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.csv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/metrics_summary.csv
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig.fasta.fai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/vdj_contig_info.pb
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/vdj_reference
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/concat_ref.bam.bai
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig.fasta
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/consensus.bam
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.bed
# /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/web_summary.html
# /Work_bio/lhn_work/cellranger/BCR/humanized_rat/refdata/reference.json
# /Work_bio/lhn_work/cellranger/BCR/humanized_rat/refdata/fasta
# /Work_bio/lhn_work/cellranger/BCR/human/refdata/reference.json
# /Work_bio/lhn_work/cellranger/BCR/human/refdata/fasta
# /Work_bio/lhn_work/cellranger/BCR/mouse/refdata/reference.json
# /Work_bio/lhn_work/cellranger/BCR/mouse/refdata/fasta
# /Work_bio/lhn_work/cellranger/BCR/humanized_mouse/refdata/reference.json
# /Work_bio/lhn_work/cellranger/BCR/humanized_mouse/refdata/fasta
# Read work_for_ath.sh
# Bash Inspect job logs, reference data sizes, and output metrics
# IN
# echo "--- test.o (tail) ---"; tail -50 /Work_bio/lhn_work/cellranger/Demo2/test.o; echo; echo "--- test.e (tail) ---"; tail -30 /Work_bio/lhn_work/cellranger/Demo2/test.e; echo; echo "--- BCR refdata sizes ---"; du -sh /Work_bio/lhn_work/cellranger/BCR/*/refdata 2>&1; echo; echo "--- metrics_summary ---"; cat /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/metrics_summary.csv

# OUT
# --- test.o (tail) ---
# 2026-06-06 06:59:00 [runtime] (join_complete)   ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.VDJ_B_GEM_WELL_PROCESSOR.LIB_VDJ_ANALYZER.SUMMARIZE_VDJ_REPORTS
# 2026-06-06 06:59:00 [runtime] (ready)           ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.VDJ_B_GEM_WELL_PROCESSOR.BUILD_PER_SAMPLE_VDJ_WS_CONTENTS
# 2026-06-06 06:59:00 [runtime] (ready)           ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.SC_MULTI_CORE.WRITE_MULTI_WEB_SUMMARY
# 2026-06-06 07:00:42 [runtime] (chunks_complete) ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.VDJ_B_GEM_WELL_PROCESSOR.LIB_VDJ_ANALYZER.VLOUPE_PREPROCESS
# 2026-06-06 07:00:42 [runtime] (ready)           ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.PICK_VDJ_OUTS
# 2026-06-06 07:00:42 [runtime] (run:local)       ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.PICK_VDJ_OUTS.fork0.chnk0.main
# 2026-06-06 07:00:43 [runtime] (chunks_complete) ID.N2-SI_TT_D1_23575YLT4.SC_VDJ_ASSEMBLER_CS.PICK_VDJ_OUTS

# Outputs:
# - Run summary HTML:                                 /work/N2-SI_TT_D1_23575YLT4/outs/web_summary.html
# - Run summary CSV:                                  /work/N2-SI_TT_D1_23575YLT4/outs/metrics_summary.csv
# - Clonotype info:                                   /work/N2-SI_TT_D1_23575YLT4/outs/clonotypes.csv
# - Filtered contig sequences FASTA:                  /work/N2-SI_TT_D1_23575YLT4/outs/filtered_contig.fasta
# - Filtered contig sequences FASTQ:                  /work/N2-SI_TT_D1_23575YLT4/outs/filtered_contig.fastq
# - Filtered contigs (CSV):                           /work/N2-SI_TT_D1_23575YLT4/outs/filtered_contig_annotations.csv
# - All-contig FASTA:                                 /work/N2-SI_TT_D1_23575YLT4/outs/all_contig.fasta
# - All-contig FASTA index:                           /work/N2-SI_TT_D1_23575YLT4/outs/all_contig.fasta.fai
# - All-contig FASTQ:                                 /work/N2-SI_TT_D1_23575YLT4/outs/all_contig.fastq
# - Read-contig alignments:                           /work/N2-SI_TT_D1_23575YLT4/outs/all_contig.bam
# - Read-contig alignment index:                      /work/N2-SI_TT_D1_23575YLT4/outs/all_contig.bam.bai
# - All contig annotations (JSON):                    /work/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.json
# - All contig annotations (BED):                     /work/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.bed
# - All contig annotations (CSV):                     /work/N2-SI_TT_D1_23575YLT4/outs/all_contig_annotations.csv
# - Barcodes that are declared to be targetted cells: /work/N2-SI_TT_D1_23575YLT4/outs/cell_barcodes.json
# - Clonotype consensus FASTA:                        /work/N2-SI_TT_D1_23575YLT4/outs/consensus.fasta
# - Clonotype consensus FASTA index:                  /work/N2-SI_TT_D1_23575YLT4/outs/consensus.fasta.fai
# - Contig-consensus alignments:                      /work/N2-SI_TT_D1_23575YLT4/outs/consensus.bam
# - Contig-consensus alignment index:                 /work/N2-SI_TT_D1_23575YLT4/outs/consensus.bam.bai
# - Clonotype consensus annotations (CSV):            /work/N2-SI_TT_D1_23575YLT4/outs/consensus_annotations.csv
# - Concatenated reference sequences:                 /work/N2-SI_TT_D1_23575YLT4/outs/concat_ref.fasta
# - Concatenated reference index:                     /work/N2-SI_TT_D1_23575YLT4/outs/concat_ref.fasta.fai
# - Contig-reference alignments:                      /work/N2-SI_TT_D1_23575YLT4/outs/concat_ref.bam
# - Contig-reference alignment index:                 /work/N2-SI_TT_D1_23575YLT4/outs/concat_ref.bam.bai
# - Loupe V(D)J Browser file:                         /work/N2-SI_TT_D1_23575YLT4/outs/vloupe.vloupe
# - V(D)J reference:
#     fasta:
#       regions: /work/N2-SI_TT_D1_23575YLT4/outs/vdj_reference/fasta/regions.fa
#     reference: /work/N2-SI_TT_D1_23575YLT4/outs/vdj_reference/reference.json
# - AIRR Rearrangement TSV:                           /work/N2-SI_TT_D1_23575YLT4/outs/airr_rearrangement.tsv
# - All contig info (ProtoBuf format):                /work/N2-SI_TT_D1_23575YLT4/outs/vdj_contig_info.pb
# - Inferred germline sequences:                      /work/N2-SI_TT_D1_23575YLT4/outs/donor_regions.fa

# Alerts:
# Not running cell annotation as there is no gene expression library in sample!

# Waiting 6 seconds for UI to do final refresh.
# Pipestance completed successfully!

# 2026-06-06 07:00:50 Shutting down.
# Saving pipestance info to "/work/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.mri.tgz"

# --- test.e (tail) ---
# Q40 bases: 5441545098(97.1527%)

# Read2 after filtering:
# total reads: 37559174
# total bases: 5601025059
# Q20 bases: 5543880993(98.9798%)
# Q30 bases: 5336151391(95.271%)
# Q40 bases: 5336151391(95.271%)

# Filtering result:
# reads passed filter: 75118348
# reads failed due to low quality: 682642
# reads failed due to too many N: 168154
# reads failed due to too short: 1749196
# reads with adapter trimmed: 4541996
# bases trimmed due to adapters: 201480656

# Duplication rate: 22.9788%

# Insert size peak (evaluated by paired-end reads): 256

# JSON report: /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.json
# HTML report: /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.html

# fastp --n_base_limit 1 --length_required 100 --in1 /fq/N2-SI_TT_D1_23575YLT4_S3_L003_R1_001.fastq.gz --out1 /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4_S2_L002_R1_001.fastq.gz --in2 /fq/N2-SI_TT_D1_23575YLT4_S3_L003_R2_001.fastq.gz --out2 /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4_S2_L002_R2_001.fastq.gz -j /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.json -h /work/01.QC/N2-SI_TT_D1_23575YLT4/N2-SI_TT_D1_23575YLT4.html -R N2-SI_TT_D1_23575YLT4 
# fastp v1.0.1, time used: 178 seconds


# tar: Removing leading `/' from member names
# tar: Removing leading `/' from hard link targets

# --- BCR refdata sizes ---
# 12M	/Work_bio/lhn_work/cellranger/BCR/humanized_mouse/refdata
# 12M	/Work_bio/lhn_work/cellranger/BCR/humanized_rat/refdata
# 12M	/Work_bio/lhn_work/cellranger/BCR/human/refdata
# 13M	/Work_bio/lhn_work/cellranger/BCR/mouse/refdata

# --- metrics_summary ---
# Estimated Number of Cells,Mean Read Pairs per Cell,Number of Cells With Productive V-J Spanning Pair,Number of Read Pairs,Valid Barcodes,Q30 Bases in Barcode,Q30 Bases in RNA Read 1,Q30 Bases in RNA Read 2,Q30 Bases in UMI,Reads Mapped to Any V(D)J Gene,Reads Mapped to IGH,Reads Mapped to IGK,Reads Mapped to IGL,Mean Used Read Pairs per Cell,Fraction Reads in Cells,Median IGH UMIs per Cell,Median IGK UMIs per Cell,Median IGL UMIs per Cell,Cells With Productive V-J Spanning Pair,"Cells With Productive V-J Spanning (IGK, IGH) Pair","Cells With Productive V-J Spanning (IGL, IGH) Pair",Paired Clonotype Diversity,Cells With IGH Contig,Cells With IGK Contig,Cells With IGL Contig,Cells With CDR3-annotated IGH Contig,Cells With CDR3-annotated IGK Contig,Cells With CDR3-annotated IGL Contig,Cells With V-J Spanning IGH Contig,Cells With V-J Spanning IGK Contig,Cells With V-J Spanning IGL Contig,Cells With Productive IGH Contig,Cells With Productive IGK Contig,Cells With Productive IGL Contig
# "2,267","16,568","1,722","37,559,174",83.0%,95.5%,97.4%,95.3%,96.6%,73.9%,22.2%,51.7%,0.0%,"3,771",49.9%,60.0,1302.0,0.0,76.0%,76.0%,0.0%,217.06,80.5%,98.7%,0.0%,79.6%,98.7%,0.0%,79.4%,98.7%,0.0%,77.4%,98.5%,0.0%
# 现在我对这个项目有了清晰的了解，给你总结一下：

# 项目概览：/Work_bio/lhn_work/cellranger
# 这是同事搭建的一套单细胞 BCR（B细胞受体）V(D)J 测序分析流水线，基于 10x Genomics 的 Cell Ranger 10.0.0 软件，封装为 Docker 容器化的可重复运行方案。

# 目录结构与功能
# 目录/文件	作用
# cellranger_10.0.0.tar	Cell Ranger 10.0.0 的 Docker 镜像离线包（约 4.6GB）
# BCR/	四个物种/品系的 V(D)J 参考数据库：human、mouse、humanized_mouse、humanized_rat（各约12-13MB）
# Demo2/	测试用的示例流程，包含输入 fastq、QC 结果、Cell Ranger 输出
# Demo2/work_for_ath.sh	核心运行脚本，定义了完整两步分析流程
# Readme	简要说明：目录介绍 + 所需输入参数 + 脚本引用位置
# 核心分析流程（work_for_ath.sh）
# Step1 - QC（质控/过滤）：用 Docker 容器内的 fastp 对 R1/R2 原始 fastq 做质控（去除含N过多的reads、长度过滤、生成 QC 报告 html/json）。

# Step2 - Cell Ranger VDJ 分析：用容器内的 cellranger vdj 命令，针对指定物种的参考库，对 QC 后的 fastq 做免疫库（IG，即 BCR 重链/轻链）的组装与注释，输出：

# 网页摘要报告 web_summary.html、指标汇总 metrics_summary.csv
# Clonotype（克隆型）信息、CDR3 注释
# contig 序列（FASTA/FASTQ/BAM）及比对结果
# Loupe V(D)J 浏览器文件 vloupe.vloupe
# AIRR 标准格式输出等
# 设计要点
# 通过 species 参数（human / mouse / humanized_mouse / humanized_rat）灵活切换参考数据库，适配人源化动物模型实验。
# 整套流程 Docker 化，避免环境依赖问题，便于在不同服务器复现。
# Demo2 中已有一个跑通的示例（样本 N2-SI_TT_D1_23575YLT4），从测试日志看：约2267个细胞，95%+碱基质量(Q30)，流程跑通成功（"Pipestance completed successfully!"）。
# 简单说：这是一套"输入测序原始数据 → 自动 QC → 调用 Cell Ranger 做单细胞 BCR 免疫库分析 → 输出克隆型/抗体序列报告"的标准化分析管线，主要用于免疫相关的单细胞测序项目（如抗体药物研发中的 B 细胞克隆分析）。
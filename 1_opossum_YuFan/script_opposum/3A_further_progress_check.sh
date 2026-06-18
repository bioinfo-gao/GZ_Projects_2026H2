# 1. 看 STAR 进程用了多少 CPU 和内存
ps aux | grep STAR | grep -v grep

# 2. 看系统 I/O wait（wa% 高说明磁盘是瓶颈）
top -bn1 | head -5

# 3. 看 STAR 用了几个线程
grep "runThreadN\|--runThread" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/f4/d48259a5e447f00d38529cc72c029f/.command.sh



tail -3 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/eb/0c3861f51392cc2188872eb71bec30/pi5_4.Log.progress.out
# Jun 17 21:04:01      4.7    69370545      295    87.4%    290.3     0.6%     3.0%     0.0%     0.0%     9.6%     0.0%_4.Log.progress.out
# Jun 17 21:05:28      4.7    69620724      295    87.4%    290.3     0.6%     3.0%     0.0%     0.0%     9.6%     0.0%
# ALL DONE!

tail -3 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/8c/bfffa0314b7b899be50ccf24bad072/NC_4.Log.progress.out
# Jun 17 21:05:28      5.8    53291668      296    87.7%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%4.Log.progress.out
# Jun 17 21:06:30      5.8    53417537      296    87.7%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%
# Jun 17 21:07:48      5.8    53543505      296    87.7%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%

# tail -3 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/8c/bfffa0314b7b899be50ccf24bad072/NC_4.Log.progress.out
# Jun 17 22:48:39      5.9    64079567      296    87.8%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%
# Jun 17 22:49:41      5.9    64205446      296    87.8%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%
# Jun 17 22:51:14      5.8    64289524      296    87.8%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%

# [eb/0c3861] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (pi5_4)                                           [ 62%] 5 of 8
# [-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT                        [  0%] 0 of 5
# [-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX                       -
# [-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS    -
# [-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT -
# [-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS -
# [-        ] NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SALMON_QUANT                                       [  0%] 0 of 5

# 状态更新
# 样本	状态	进度
# pi5_4	✅ 完成	ALL DONE!
# NC_4	🔄 运行中	64.3M / ~70M reads (92%)
# NC_4 预计还需 约 1 小时完成（5.8M reads/hr，剩余 ~5.7M reads）。

# NC_4 完成后，最后一批 STAR（pi5_1 + NC_1）会自动启动，预计再需 8-12 小时。

# STAR_ALIGN 总进度：7/8 即将完成。

# tail -3 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/8c/bfffa0314b7b899be50ccf24bad072/NC_4.Log.progress.out
# Jun 18 00:12:01      5.9    72765918      296    87.8%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%4.Log.progress.out
# Jun 18 00:13:27      5.9    72933831      296    87.8%    290.6     0.7%     3.2%     0.0%     0.0%     9.0%     0.0%
# ALL DONE!

-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX                       -
executor >  local (8)
[54/b347ff] NFCORE_RNASEQ:PREPARE_GENOME:GTF_FILTER (dv-2k.fasta)                                        [100%] 1 of 1, cached: 1 ✔                                                                                                      "us1" 23:46 17-Jun-26
[05/039a6b] NFCORE_RNASEQ:PREPARE_GENOME:GTF2BED (dv-2k.filtered.gtf)                                    [100%] 1 of 1, cached: 1 ✔
[c5/ba0a82] NFCORE_RNASEQ:PREPARE_GENOME:MAKE_TRANSCRIPTS_FASTA (rsem/dv-2k.fasta)                       [100%] 1 of 1, cached: 1 ✔
[c9/ca4b36] NFCORE_RNASEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES (dv-2k.fasta)                              [100%] 1 of 1, cached: 1 ✔
[-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:CAT_FASTQ                          -
[55/2194e6] NFC…SEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC (pi5_1) [100%] 8 of 8, cached: 8 ✔
[a5/b4b52f] NFC…FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE (pi5_4) [100%] 8 of 8, cached: 8 ✔
[-        ] NFC…NASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_INDEX -
[-        ] NFC…NASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_SUBSAMPLE_FQ_SALMON:FQ_SUBSAMPLE -
[-        ] NFC…NASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_QUANT -
[8c/bfffa0] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (NC_4)                                            [ 75%] 6 of 8
[-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT                        [  0%] 0 of 6
[-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX                       -
[-        ] NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS    -
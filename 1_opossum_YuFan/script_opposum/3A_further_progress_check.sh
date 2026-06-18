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


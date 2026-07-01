# zcat your_R1.fastq.gz | head -n 40 | sed -n '2~4p'
zcat /home/gao/Dropbox/Quote_06062601_output/OVO/OVO_CKDL260011462-1A_23JCTGLT3_L7_1.fq.gz | head -n 4 
zcat /home/gao/Dropbox/Quote_06062601_output/OVO/OVO_CKDL260011462-1A_23JCTGLT3_L7_1.fq.gz | head -n 40 | sed -n '2~4p'
zcat /home/gao/Dropbox/Quote_06062601_output/OVO/OVO_CKDL260011462-1A_23JCTGLT3_L7_2.fq.gz | head -n 40 | sed -n '2~4p'

1. 模板匹配分析
你的数据中包含了不同的 Template 类型，我们可以清晰地通过 Read 1 的开头识别出来：
Read 3, 4, 7, 8 匹配 Template 1 (3N + AGTCGA)：
以 Read 3 为例：GNC (3bp N) + AGTCGA (固定) + CCTATTTTCAATTTAACGTCG (骨架) + GCAATGAACAAGAGCAGCGA (sgRNA) + GTTTTAG...
结构完美吻合。
Read 2, 6, 9 匹配 Template 2 (4N + TCAGCT)：
以 Read 2 为例：GNGA (4bp N) + TCAGCT (固定) + CCTATTTTCAATTTAACGTCG (骨架) + CCCGTGCTGGTCATGTCCAG (sgRNA) + GTTTTAG...
结构完美吻合。
Read 1, 5, 10 是背景噪音：
这些序列看起来像基因组污染或者是没扩增好的副产物，这在 CRISPR 库测序中是正常的，通常占 5-15% 左右。
2. 为什么不需要读 R2？
从你给出的 R1 来看，这一条序列里已经包含了：
UMI (前面的 N)： 用于去重。
Sample ID (固定碱基)： 用于区分是 Template 1、2 还是 3。
sgRNA (20bp)： 位于骨架序列 ...ACGTCG 之后。这是你筛选实验的核心数据。
R2 通常是 R1 的反向互补，在扩增子测序中，R2 往往只能测到 sgRNA 的末端或者是 Scaffold 区域，提供的信息通常是冗余的。分析这个项目，R1 是绝对的主角。







zcat /home/gao/Dropbox/Quote_06062601_output/Undetermined/Undetermined_Undetermined_23JCTGLT3_L7_1.fq.gz | head -n 40 | sed -n '2~4p'
zcat /home/gao/Dropbox/Quote_06062601_output/Undetermined/Undetermined_Undetermined_23JCTGLT3_L7_2.fq.gz | head -n 40 | sed -n '2~4p'
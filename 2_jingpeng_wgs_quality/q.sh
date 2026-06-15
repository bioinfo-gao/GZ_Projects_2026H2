/home/gao/Dropbox/JinPeng]:
gao@us1 $ tree -l
.
├── 02.Report_X202SC25091740-Z01-F012.zip
├── 23TF7FLT4_3_0469165296_8464_S1_L003_R1_001.fastq.gz
├── 23TF7FLT4_3_0469165296_8464_S1_L003_R2_001.fastq.gz
├── 23TF7FLT4_3_0469165296_8465_S2_L003_R1_001.fastq.gz
├── 23TF7FLT4_3_0469165296_8465_S2_L003_R2_001.fastq.gz
├── 23TF7FLT4_5_0469165288_38681st_S9_L005_R1_001.fastq.gz
├── 23TF7FLT4_5_0469165288_38681st_S9_L005_R2_001.fastq.gz
├── checkSize.xls
├── J_3685
│   ├── J_3685_CKDL260010985-1A_23JJWKLT3_L3_1.fq.gz
│   ├── J_3685_CKDL260010985-1A_23JJWKLT3_L3_2.fq.gz
│   └── MD5.txt
├── J_3689
│   ├── J_3689_CKDL260010985-1A_23JJWKLT3_L3_1.fq.gz
│   ├── J_3689_CKDL260010985-1A_23JJWKLT3_L3_2.fq.gz
│   └── MD5.txt
├── J_3852R0B1TP
│   ├── J_3852R0B1TP_CKDL260010985-1A_23JJWKLT3_L3_1.fq.gz
│   ├── J_3852R0B1TP_CKDL260010985-1A_23JJWKLT3_L3_2.fq.gz
│   └── MD5.txt
├── J_3852R0B2TP
│   ├── J_3852R0B2TP_CKDL260010985-1A_23JJWKLT3_L3_1.fq.gz
│   ├── J_3852R0B2TP_CKDL260010985-1A_23JJWKLT3_L3_2.fq.gz
│   └── MD5.txt
├── J_3852R0origin
│   ├── J_3852R0origin_CKDL260010985-1A_23JJWKLT3_L3_1.fq.gz
│   ├── J_3852R0origin_CKDL260010985-1A_23JJWKLT3_L3_2.fq.gz

# kimi

# 我来帮你设计一个快速检查方案，针对酶切时间不够可能导致的样品质量问题。
# 核心问题：酶切时间不够的影响
# 表格
# 问题	           表现
# 酶切不完全	基因组 DNA 残留多，RNA-seq 中 rRNA/基因组比例升高
# RNA降解       片段大小变小，RIN 值低
# 接头二聚体	小片段 (<100bp) 比例升高
# 复杂度降低	重复序列多，唯一比对率下降


# 快速抽样检查方案（推荐）
bash
#!/bin/bash
# fastq_quality_check.sh - 快速检查酶切样品质量

DIR="/home/gao/Dropbox/JinPeng"
SAMPLES=$(ls -d $DIR/J_* 2>/dev/null | head -5)  # 先抽5个
OUTPUT="qc_report_$(date +%Y%m%d).txt"

echo "=== FastQ 质量抽样检查 ===" > $OUTPUT
echo "时间: $(date)" >> $OUTPUT
echo "" >> $OUTPUT

for sample_dir in $SAMPLES; do
    sample=$(basename $sample_dir)
    r1="${sample_dir}/${sample}_*_L3_1.fq.gz"
    r2="${sample_dir}/${sample}_*_L3_2.fq.gz"
    
    echo "=== $sample ===" >> $OUTPUT
    
    # 1. 文件大小检查（异常小可能降解严重）
    size1=$(ls -lh $r1 2>/dev/null | awk '{print $5}')
    size2=$(ls -lh $r2 2>/dev/null | awk '{print $5}')
    echo "R1大小: $size1, R2大小: $size2" >> $OUTPUT
    
    # 2. 读取总条数（异常少可能有问题）
    count1=$(zcat $r1 2>/dev/null | wc -l)
    reads1=$((count1 / 4))
    echo "总reads: $reads1" >> $OUTPUT
    
    # 3. 抽样检查序列长度分布（关键！）
    echo "序列长度分布 (前1000条):" >> $OUTPUT
    zcat $r1 2>/dev/null | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -k2 -n | tail -5 >> $OUTPUT
    
    # 4. 检查短片段比例（<<50bp 可能是降解或二聚体）
    short=$(zcat $r1 2>/dev/null | awk 'NR%4==2 {if(length($0)<50) count++} END{print count+0}')
    total=$(zcat $r1 2>/dev/null | awk 'NR%4==2 {count++} END{print count}')
    short_pct=$(awk "BEGIN{printf \"%.2f\", $short/$total*100}")
    echo "短片段(<50bp)比例: $short_pct%" >> $OUTPUT
    
    # 5. 简单碱基质量检查（Q30比例）
    echo "Q30比例 (抽样):" >> $OUTPUT
    zcat $r1 2>/dev/null | head -4000 | awk 'NR%4==0 {
        q30=0; total=0;
        for(i=1;i<=length($0);i++) {
            q=substr($0,i,1);
            if(q!="") {
                total++;
                if((q-33)>=30) q30++;
            }
        }
        if(total>0) print q30/total*100
    }' | awk '{sum+=$1; count++} END{if(count>0) printf "Q30: %.1f%%\n", sum/count}' >> $OUTPUT
    
    echo "" >> $OUTPUT
done

echo "报告已保存: $OUTPUT"
cat $OUTPUT

# 更快的方法：用现有工具
# bash
# # 方法1: seqkit 快速统计（如果已安装）
# for d in J_*; do
#     echo "=== $d ==="
#     seqkit stats ${d}/*.fq.gz 2>/dev/null | head -2
# done

# # 方法2: 检查MD5确保文件完整
# for d in J_*; do
#     (cd $d && md5sum -c MD5.txt 2>/dev/null) && echo "$d: MD5 OK" || echo "$d: MD5 FAIL"
# done

# # 方法3: 快速看reads数和文件大小
# du -sh J_*/*
# 关键判断指标
# 表格
# 指标	正常范围	异常提示
# 文件大小	与同类样品一致	明显偏小（<<50%）
# Reads数	>10M	<5M 需警惕
# 短片段<<50bp	<5%	>20% 可能酶切过度或降解
# Q30比例	>85%	<75% 质量差
# 序列长度	稳定150bp	长度分布宽，很多<<100bp
# 我的建议
# 抽样检查就足够了，原因：
# 同一批次酶切条件相同，质量问题应系统性出现
# 抽5个样品能看出趋势
# 真正关键的是后续比对后的指标（rRNA残留率、基因组比对率、重复率）
# 抽样后如果发现问题，再全部检查；如果正常，直接跑流程，看比对结果。
# 你想运行哪个方案？或者先检查 J_3685 和 J_3689 两个看看？
# 高峰时段算力不足，已切换至 K2.6 快速，升级会员畅用思考模型

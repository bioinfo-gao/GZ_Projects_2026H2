# 1. 查看实际运行的样本列表
cat /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/nf_core_samplesheet.csv

# 2. 统计 STAR_ALIGN 任务总数
grep "STAR_ALIGN" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/.nextflow.log | grep -c "Submitted process"

# 3. 查看已完成和运行中的任务
grep "STAR_ALIGN" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/.nextflow.log | grep -E "COMPLETED|RUNNING" | wc -l
# 1664

# 1. 统计 samplesheet 中的样本数
tail -n +2 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/nf_core_samplesheet.csv | wc -l
# 预期输出: 34

# 2. 查看哪些样本的 STAR_ALIGN 已完成
grep "STAR_ALIGN" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/.nextflow.log | \
  grep "COMPLETED" | \
  grep -oP "STAR_ALIGN \(\K[^)]+" | \
  sort -u | wc -l
# 预期输出: ~33

# 3. 查看哪个样本还在运行
grep "STAR_ALIGN" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/.nextflow.log | \
  grep "RUNNING" | \
  grep -oP "STAR_ALIGN \(\K[^)]+"
# 预期输出: 剩余 1 个样本名

# OMPLETED # 16
# 本次启动后真正运行并成功完成的任务
# 您刚运行的 grep "COMPLETED" 结果
# Cached # 17
# 之前运行过，Nextflow 直接复用结果的任务
# 日志中大量的 Cached process > ... 行
# Total Done # 33
# 16 + 17 = 33
# 进度条显示的 33 of 34
# RUNNING # 1
# 剩下的最后一个任务 (Wh2b)
# 进度条显示 97%

# 剩余 1 个任务正在运行中（即 Wh2b）
# # 1. 定位 Wh2b 最新的工作目录
# WH2B_DIR=$(find /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work -maxdepth 3 -name ".command.run" -exec grep -l "Wh2b" {} + | xargs ls -td | head -n 1 | xargs dirname)

# # 2. 查看 STAR 进度文件（通常位于工作目录根目录或 STAR/ 子目录）
# tail -n 15 ${WH2B_DIR}/Log.progress.out 2>/dev/null || tail -n 15 ${WH2B_DIR}/STAR/Log.progress.out

# 更高效的进度查看方式
# 1. Nextflow 原生进度查看（推荐）
# 查看当前运行会话的实时任务状态（无需猜目录）
#nextflow log -name compassionate_nightingale -f name,status,workdir,complete | grep Wh2b
# cd /home/gao/projects/2026_Item9_gc/output_results/star_salmon/log 
# cd /home/gao/projects/2026_Item9_gc/scripts 
# # It looks like no pipeline was executed in this folder (or execution history is empty)
# nextflow log compassionate_nightingale -f name,status,workdir,complete | grep Wh2b # It looks like no pipeline was executed in this folder (or execution history is empty)



# 在项目目录下执行

grep "Wh2b" /home/gao/projects/2026_Item9_gc/output_results/pipeline_info/execution_trace_*.txt >> wh2b_progress.txt

# 2. 批量查看所有样本 STAR 进度
for sample in WTS4_1 WTS4_2 WTS4_4 WTS4_5 WhS4_1 WhS4_4 WhS4_5 WhS4_9 Wh1a Wh1b Wh1c Wh1d Wh2a Wh2b Wh2c Wh2d Wh3a Wh3b Wh3c Wh3d Wh4a Wh4b Wh4c Wh4d Wh5a Wh5b Wh5c Wh5d Wh6a Wh6b Wh6c Wh6d TperMix TtriMix; do
  dir=$(find /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work -name ".command.run" -exec grep -l "$sample" {} + 2>/dev/null | head -1 | xargs dirname)
  if [ -f "$dir/Log.progress.out" ]; then
    pct=$(tail -1 "$dir/Log.progress.out" | grep -oP '% mapped \| \K[0-9.]+')
    echo -e "$sample\tMapped: ${pct}%"
  fi
done | sort -k2 -n

# Wh2b 已经是那 33 个已完成的样本之一！
# 查看当前运行中的 STAR_ALIGN 任务
grep "STAR_ALIGN.*RUNNING" /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/.nextflow.log | tail -5

# 或查看执行追踪中未完成的任务
grep "STAR_ALIGN" /home/gao/projects/2026_Item9_gc/output_results/pipeline_info/execution_trace_*.txt | grep -v "COMPLETED\|CACHED" | tail -10



# 找到 WTS4_1 的工作目录
WTS4_1_DIR=$(find /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work -name ".command.run" -exec grep -l "WTS4_1" {} + | xargs ls -td | head -1 | xargs dirname)
echo $WTS4_1_DIR

# 查看 STAR 进度（如果还在运行）
tail -n 20 ${WTS4_1_DIR}/STAR/Log.progress.out 2>/dev/null || echo "STAR 已完成或未生成进度文件"

# 查看任务日志
tail -n 50 ${WTS4_1_DIR}/.command.log 2>/dev/null | grep -E "ERROR|WARNING|Completed|mapped"

tail -n 3 /home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/b0/b6ed769b0abaddef055149e3eaadca/Log.progress.out


# WTS4_1 的 STAR 工作目录（根据您之前日志）
WTS4_1_STAR_DIR="/home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/work/b0/b6ed769b0abaddef055149e3eaadca/"

# 查看 STAR 实时比对进度
echo "=== STAR Log.progress.out ==="
tail -n 10 ${WTS4_1_STAR_DIR}/WTS4_1.Log.progress.out 2>/dev/null || echo "⚠️ 文件不存在或未生成"


时间范围：18:50:48 → 19:02:27（12 分钟监控窗口）

列含义（从左到右）：
tail -n 10 ${NEW_WTS4_DIR}/WTS4_1.Log.progress.out
           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped
                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other
Apr 07 22:38:36 Started 1st pass mapping

Apr 08 18:50:48      4.1    10019367      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:51:48      4.1    10061470      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:53:17      4.1    10145592      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:54:21      4.1    10271866      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:55:24      4.1    10313943      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:57:00      4.1    10398126      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.2%
Apr 08 18:58:00      4.1    10524363      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 18:59:08      4.1    10608495      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 19:00:38      4.1    10650562      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 19:02:27      4.1    10818863      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%

tail -n 10 ${WTS4_1_STAR_DIR}/WTS4_1.Log.progress.out 2>/dev/null || echo "⚠️ 文件不存在或未生成"
Apr 08 19:10:08      4.1    11324241      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 19:12:10      4.1    11408820      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 19:13:10      4.1    11535667      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.5%     4.1%
Apr 08 19:14:25      4.1    11577912      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.6%     4.1%
Apr 08 19:16:17      4.1    11662442      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.6%     4.1%
Apr 08 19:18:20      4.1    11831406      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.6%     4.1%# 查看当前运行的进程
ps aux | grep samtools

# 查看内存和 CPU 使用
top -p $(pgrep samtools)

# 或者使用 htop 更直观地查看
htop
Apr 08 19:20:19      4.1    11915884      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.6%     4.2%
Apr 08 19:22:28      4.1    12084778      294    71.9%    294.9     0.3%    10.2%    12.2%     0.0%     1.6%     4.2%
Apr 08 19:24:28      4.0    12169222      294    71.9%    294.8     0.3%    10.2%    12.2%     0.0%     1.6%     4.2%
Apr 08 19:26:25      4.1    12338081      294    71.9%    294.8     0.3%    10.2%    12.2%     0.0%     1.6%     4.2%


# 如果文件存在，解析 % mapped 字段
if [ -f "${WTS4_1_STAR_DIR}/Log.progress.out" ]; then
    echo ""
    echo "=== 比对完成比例 ==="
    tail -1 ${WTS4_1_STAR_DIR}/Log.progress.out | grep -oP '% mapped \| \K[0-9.]+' | xargs -I {} echo "已比对: {}%"
fi


# 估算 WTS4_1 的总 reads 数
zcat /hozcat /home/gao/Dropbox/Quote_03062602_plant/WTS4_1/WTS4_1_R1.fq.gz | echo $((`wc -l`/4))
22,286,484


当前进度: ~48.5%
已用时间: ~7.5 小时 (11:46 → 19:02)
处理速度: ~4.1 M/hr

剩余 reads: 22.3M - 10.8M = ~11.5M
剩余时间: 11.5M ÷ 4.1M/hr ≈ 2.8 小时

预计完成:
├── STAR 比对完成: ~21:45-22:00
├── BAM 排序/索引: ~22:15
├── 下游任务触发: ~22:30
└── 总剩余时间: **约 3-4 小时**



# 查看当前运行的进程
ps aux | grep samtools

# 查看内存和 CPU 使用
top -p $(pgrep samtools)

# 或者使用 htop 更直观地查看
htop


# 查看可用内存
free -h

# 查看磁盘空间
df -h /home/gao/work

# 查看 Nextflow 是否还在运行
ps aux | grep nextflow | grep -v grep

# 查看当前的 pipeline 状态
nextflow log -l 1


# 如果 Nextflow 仍在运行
# 那就继续等待，不要中断它。从你的输出来看，系统正在处理其他任务（RSeqC 等），SAMTOOLS_SORT 那个失败的任务可能已经被跳过了。

# ⚠️ 如果你想强制重新运行那个失败的任务
# 可以只针对特定任务重新运行，而不影响其他结果

# 找到失败任务的 hash
nextflow log -f process,task_id,status | grep FAILED

# 然后只恢复这个任务（不会影响其他已完成的工作）
nextflow resume -resume


# 查看该任务的日志
find /home/gao/projects/2026_Item9_gc/scripts/work -path "*d0c685*" -name ".command.log" -exec tail -20 {} \;
find /home/gao/projects/2026_Item9_gc/scripts/work -path "*d0c685*" -name ".command.err" -exec cat {} \;

#1:10am finished
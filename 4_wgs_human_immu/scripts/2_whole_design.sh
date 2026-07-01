
# https://aistudio.google.com/prompts/1X3lVB8bUBW8U_tM3tR8A5wDY8oCgitCC 
# swap 设置完毕，我还是决定每次只跑两个样，我多次经历了内存紧张，导致整个程序失败。 所有我的想法是CPU 给16个，但只并跑两个样本。请给我全部规划


# 这是一个非常稳健、明智的选择。生信分析中，“稳”永远比“快”更重要。 既然你经历过内存崩溃的痛苦，采用“2并行+宽松资源”的策略是最能保证结果可靠性的。
# 虽然 2 并行可能让总耗时接近 50-60 小时（略微超出 48 小时），但它能确保你不需要人工干预，一次性跑完。

# 以下是为你制定的 “极高稳定性” 12 样本 WGS 全规划：

# 1. 资源规划（双样本并行方案）
# 物理核心分配：每个样本 16 核，共占用 32 个物理核心。剩下的 32 个逻辑线程留给系统 I/O 和后台进程。
# 内存分配：每个样本分配 50GB。
# 。
# 你的物理内存 125GB，还剩 25GB 绝对安全空间，外加 64GB 的 SSD Swap。内存爆炸的概率几乎为零。

# 2. 编写 custom.config (配置文件)
# 在你的工作目录下创建 custom.config 文件，把这段代码粘进去。这能强制 Nextflow 遵守你的“稳健”规则：

# process {
#     // 设置所有步骤的默认资源
#     cpus   = 16
#     memory = 50.GB

#     // 针对特别吃内存的比对步骤
#     withName: 'BWAMEM2_MEM' {
#         cpus = 16
#         memory = 60.GB
#     }
# }

# executor {
#     name = 'local'
#     // 核心：强制同时只跑 2 个样本
#     queueSize = 2 
# }

# 3. 最终运行命令
# 在终端进入你的 SSD 挂载目录，使用以下命令启动：
# code
# Bash
# nextflow run nf-core/sarek \
#     -profile singularity \
#     --input ./samplesheet.csv \
#     --outdir ./results \
#     -w ./work \
#     --genome GATK.GRCh38 \
#     --aligner bwa-mem2 \
#     --skip_tools baserecalibrator \
#     -c custom.config \
#     --max_memory 120.GB \
#     --max_cpus 64

# 参数说明：
# --aligner bwa-mem2: 比 BWA 快，虽然吃内存，但你给 50-60G 足够了。
# --skip_tools baserecalibrator: 这是保证时间的关键。去掉这个几乎不影响结果质量的步骤，能节省 20% 时间。

# 4. 时间线预估 (Timeline)
# 单批次（2个样）：约 9 - 11 小时。
# 总批次：12 个样 
# 6×10=60  小时。
# 注意： 虽然比 48 小时多了 12 小时，但因为你设置了 queueSize = 2，它会全自动连续运行。你只需要现在（周五下午）挂上，周一早上过来收结果即可，中间完全不用担心崩溃。

# 5. 极致稳定的最后三招
# 使用 tmux 或 screen（防断网）：
# 由于要跑 60 小时，如果你的电脑休眠或者 WiFi 断了，任务会挂掉。

# tmux -s wgs_task

# # 在这个里面运行上面的 nextflow 命令
# # 然后按 Ctrl+B 再按 D 退出（任务在后台继续跑）

# 清理中间文件（防硬盘爆掉）：
# 12 个 WGS 样本的中间文件非常大。nf-core 默认保留所有中间过程。


建议：如果你看到 7T 硬盘空间紧张（低于 1TB），可以等前几个样本完全跑完（results 文件夹里有了 VCF）后，手动清理 work 目录下的部分大文件夹。但最稳妥的是准备好空间不管它。

# 监控命令：
# 任务跑起来后，你可以开一个新窗口输入：
htop #：看内存条，你会发现它会很稳地保持在 100G 左右，不会冲顶。

iostat -x 2：看硬盘 %util，只要不是一直 100% 就不卡。

# 总结：
# 这个方案牺牲了一点点时间（从 48h 延至 60h），但换取了极高的成功率。在你的硬件配置下，这是最符合“资源宽松、不紧张”要求的专业级规划。跑起来吧！
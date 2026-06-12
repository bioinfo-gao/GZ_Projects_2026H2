# 一个经典的生物信息学挑战：拥有高质量的基因组（DNA Zoo），但缺乏同物种的注释文件（GTF/GFF）。
# 对于北美负鼠（Didelphis virginiana）这种非模式生物，最好的策略是利用近缘物种（南美灰短尾负鼠 Monodelphis domestica）的注释进行“同源迁移”，并结合你手头的 RNA-seq 数据进行“自动校正”。

# 以下是推荐的三步走策略：
# 策略一：使用 Liftoff 进行注释迁移（最快、最推荐）
# Liftoff 是目前最流行的工具，它可以将近缘物种的 GFF/GTF 注释直接“投射”到你的新基因组上。
# 准备文件：
# 靶标基因组：你刚下载的 DNA Zoo D. virginiana 基因组（Fasta）。
# 参考基因组：南美灰负鼠的基因组（Fasta）。
# 参考注释：南美灰负鼠的注释文件（GFF/GTF）。
# 操作步骤：
# code
# Bash
# # 安装 Liftoff (建议在 conda 环境)
# pip install liftoff

# # 运行迁移
# liftoff -g Monodelphis_domestica.gff \
#         -o Didelphis_virginiana_target.gff \
#         Didelphis_virginiana_genome.fasta \
#         Monodelphis_domestica_genome.fasta
# 优点： 能够保留基因名称和功能信息，方便后续分析。
# 策略二：基于 RNA-seq 数据的 StringTie 自组装（最准确、反映真实表达）
# 既然你已经有了 RNA-seq 数据，为什么不让数据说话？你可以利用这些 reads 直接在 DNA Zoo 基因组上组装出属于北美负鼠的转录本。
# 标准流程（New Tuxedo Pipeline）：
# 比对： 使用 HISAT2 或 STAR 将 RNA-seq 的 reads 比对到 DNA Zoo 的基因组上（输出 BAM 文件）。
# 组装： 使用 StringTie 对每个样本进行组装。
# code
# Bash
# stringtie sample.bam -o sample.gtf -p 8
# 合并： 使用 stringtie --merge 将所有样本的 GTF 合并成一个完整、统一的北美负鼠注释文件。
# 优点： 这是最贴合你实验样本的注释，能发现北美负鼠特有的、而在南美负鼠中没有的异构体（Isoforms）或新基因。
# 策略三：混合策略（工业级标准）
# 这是目前发表高分文章最稳妥的做法：
# 先用 Liftoff 将南美负鼠的注释迁移过来作为“骨架”。
# 再用 StringTie 结合你的 RNA-seq 数据进行组装。
# 最后用 gffcompare 将两者合并。
# 已知基因沿用南美负鼠的 ID。
# 新发现的转录本标记为 novel。
# 针对你的差异分析（DE Analysis）的最终工作流建议：
# 建立索引：
# 使用 hisat2-build 为 DNA Zoo 基因组建立索引。
# 生成注释：
# 运行一次 Liftoff，把南美灰负鼠的基因映射到 DNA Zoo 基因组上，得到 Dv_ref.gtf。
# 定量：
# 使用 FeatureCounts 或 StringTie 进行表达定量。
# 如果你用了 Liftoff 的 GTF，定量后的基因 ID 就会是南美负鼠的基因名（如 TP53, GAPDH），这极其利于你后续做 GO/KEGG 功能富集分析。
# 差异分析：
# 将 count 矩阵导入 DESeq2 或 edgeR。
# 💡 管理提示
# Vim 操作： 在编写比对和定量的 Shell 脚本时，记得使用 vim run_pipeline.sh 进行编辑。
# 结果交付： 差异分析完成后产生的 PDF 图片和 CSV 报表，可以直接移动到 /srv/external_data/client_sftp/files 目录下。
# 客户查看： 你的客户可以通过 bore.pub:54833 (或者当前的实时端口) 立即通过 WinSCP 下载结果。
# 你需要我为你提供具体的 HISAT2 或 StringTie 的命令行参数代码吗？

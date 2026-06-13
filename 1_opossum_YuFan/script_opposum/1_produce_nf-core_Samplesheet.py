# 运行方法
# cd /home/gao/projects_2026H1/2026_Item7_LJZ/scripts
# python 1_produce_nf-core_Samplesheet.py
# mamba activate regular_bioinfo
import pandas as pd
import glob
import os

# 配置路径
fastq_base_dir = "/home/gao/Dropbox/0602/01.RawData/"

original_csv = "/home/gao/projects_2026H1/2026_Item16_ZhenYan/scripts/zy.csv"
output_samplesheet = "/home/gao/projects_2026H1/2026_Item16_ZhenYan/scripts/nf_core_samplesheet.csv"

# 读取原始表格
df = pd.read_csv(original_csv)
samples = df['Name in File'].unique() # 获取所有样本名称, ATH 使用这个固定格式

def convert_sample_name(csv_name):
    """Convert CSV sample name to directory name format by replacing hyphens with underscores"""
    # Replace hyphens with underscores to match directory names
    dir_name = csv_name.replace("-", "_")
    return dir_name

data = []
for sample in samples:
    # 转换样本名称以匹配目录名称
    dir_name = convert_sample_name(sample)
    # 构建样本目录路径
    sample_dir = os.path.join(fastq_base_dir, dir_name)
    
    if not os.path.exists(sample_dir):
        print(f"Warning: 找不到样本 {sample} (目录名: {dir_name}) 的目录 {sample_dir}")
        continue
    
    # 查找 R1 和 R2 文件
    # 格式: {dir_name}_CKDL..._1.fq.gz / {dir_name}_CKDL..._2.fq.gz
    r1_files = glob.glob(os.path.join(sample_dir, f"{dir_name}_*_1.fq.gz"))
    r2_files = glob.glob(os.path.join(sample_dir, f"{dir_name}_*_2.fq.gz"))
    
    r1 = None
    r2 = None
    
    # 使用找到的文件
    if r1_files and r2_files:
        r1 = r1_files[0]
        r2 = r2_files[0]
    
    # 检查文件是否存在
    if r1 and r2:
        # strandedness 设为 reverse (如果是标准库)
        # strandedness 设为 auto 让流程自动检测，或者设为 reverse (如果是标准库)
        data.append([sample, r1, r2, 'reverse'])
    else:
        print(f"Warning: 找不到样本 {sample} (目录名: {dir_name}) 的 fastq 文件于 {sample_dir}")

# 生成新表格
final_df = pd.DataFrame(data, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
final_df.to_csv(output_samplesheet, index=False)
print(f"成功生成 nf-core 专用 Samplesheet: {output_samplesheet}")
print(f"共处理 {len(data)} 个样本")

# 运行方法
# cd /home/gao/projects/2026_Item7_LJZ/scripts
# python 1_produce_nf-core_Samplesheet.py
# 必要注释，在 1_produce_nf-core_Samplesheet.py (line 22) 里主要解释了 3 件事：

# 为什么不能只取第一对 R1/R2
# 为什么要用去掉末尾 _1/_2 的共同前缀来配对
# 为什么同一个 sample 需要在 nf-core samplesheet 里输出多行

import pandas as pd
import glob
import os

# 配置路径
# fastq_base_dir = "/home/gao/Dropbox/Quote_03032601/"
fastq_base_dir = "/home/gao/Dropbox/0505/01.RawData/rRNADepletionKitTest"

original_csv = "/home/gao/projects/2026_Item12_rRNA/scripts/samples.csv"
output_samplesheet = "/home/gao/projects/2026_Item12_rRNA/scripts/nf_core_samplesheet.csv"

def convert_sample_name(csv_name):
    """Convert CSV sample name to directory name format by replacing hyphens with underscores"""
    # Replace hyphens with underscores to match directory names
    dir_name = csv_name.replace("-", "_")
    return dir_name

def build_fastq_pairs(sample_dir, dir_name):
    """Find all paired FASTQ files for one sample directory."""
    # 同一样品可能被拆成多次测序，因此这里不能只取第一对 R1/R2。
    r1_files = sorted(glob.glob(os.path.join(sample_dir, f"{dir_name}_*_1.fq.gz")))
    r2_files = sorted(glob.glob(os.path.join(sample_dir, f"{dir_name}_*_2.fq.gz")))

    # 用去掉末尾 _1 / _2 的共同前缀做 key，确保只配对同一批次的 reads。
    r1_map = {
        os.path.basename(path).replace("_1.fq.gz", ""): path
        for path in r1_files
    }
    r2_map = {
        os.path.basename(path).replace("_2.fq.gz", ""): path
        for path in r2_files
    }

    shared_keys = sorted(set(r1_map) & set(r2_map))
    missing_r2 = sorted(set(r1_map) - set(r2_map))
    missing_r1 = sorted(set(r2_map) - set(r1_map))

    for key in missing_r2:
        print(f"Warning: {sample_dir} 中 {key} 缺少对应的 R2 文件")
    for key in missing_r1:
        print(f"Warning: {sample_dir} 中 {key} 缺少对应的 R1 文件")

    # nf-core samplesheet 支持同一个 sample 写多行，每行对应一对 FASTQ。
    return [(r1_map[key], r2_map[key]) for key in shared_keys]

def main():
    # 读取原始表格
    df = pd.read_csv(original_csv)
    samples = df['Name in File'].unique() # 获取所有样本名称, ATH 使用这个固定格式

    data = []
    for sample in samples:
        # 转换样本名称以匹配目录名称
        dir_name = convert_sample_name(sample)
        # 构建样本目录路径
        sample_dir = os.path.join(fastq_base_dir, dir_name)

        if not os.path.exists(sample_dir):
            print(f"Warning: 找不到样本 {sample} (目录名: {dir_name}) 的目录 {sample_dir}")
            continue

        fastq_pairs = build_fastq_pairs(sample_dir, dir_name)

        if fastq_pairs:
            for r1, r2 in fastq_pairs:
                # 同一样品若有多对 FASTQ，这里会追加多行，而不是覆盖前一对。
                # strandedness 设为 reverse (如果是标准库)
                # strandedness 设为 auto 让流程自动检测，或者设为 reverse (如果是标准库)
                data.append([sample, r1, r2, 'reverse'])
        else:
            print(f"Warning: 找不到样本 {sample} (目录名: {dir_name}) 的 fastq 文件于 {sample_dir}")

    # 生成新表格
    final_df = pd.DataFrame(data, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
    final_df.to_csv(output_samplesheet, index=False)
    print(f"成功生成 nf-core 专用 Samplesheet: {output_samplesheet}")
    print(f"共生成 {len(data)} 行 paired-end 记录")

if __name__ == "__main__":
    main()

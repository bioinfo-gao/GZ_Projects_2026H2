#!/usr/bin/env python3
# 运行方法：cd /home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts && conda run -n regular_bioinfo python 1_produce_nf-core_Samplesheet.py
#
# 注意: 该项目 FASTQ 文件名使用 "Name in File" 列的值 (如 cDC1_TNFa_1)，
# 而不是 Group+rep 组合 (Group="TNFa" 在文件名中实际是 "cDC1_TNFa")。
import pandas as pd
import glob
import os

fastq_base_dir = "/Work_bio/dropbox/Dropbox_Data/Quote_06032601_Lijian_Wu/"
original_sheet = "../Sample_Sheet_Lijian_Wu.xlsx"
output_samplesheet = "nf_core_samplesheet.csv"

df = pd.read_excel(original_sheet)
df = df[df['Group'].notna()].reset_index(drop=True)

data = []
for _, row in df.iterrows():
    name_in_file = str(row['Name in File']).strip()
    sample_id = name_in_file

    pattern_r1 = os.path.join(fastq_base_dir, f"*_{name_in_file}_S*_L*_R1_001.fastq.gz")
    pattern_r2 = os.path.join(fastq_base_dir, f"*_{name_in_file}_S*_L*_R2_001.fastq.gz")
    r1_files = sorted(glob.glob(pattern_r1))
    r2_files = sorted(glob.glob(pattern_r2))

    if r1_files and r2_files:
        data.append([sample_id, r1_files[0], r2_files[0], 'reverse'])
    else:
        print(f"WARNING: FASTQ files not found for {sample_id} (pattern: {pattern_r1})")

final_df = pd.DataFrame(data, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
final_df.to_csv(output_samplesheet, index=False)
print(f"Samplesheet saved: {output_samplesheet}")
print(f"Total samples: {len(data)} / {len(df)}")
print(final_df.to_string())

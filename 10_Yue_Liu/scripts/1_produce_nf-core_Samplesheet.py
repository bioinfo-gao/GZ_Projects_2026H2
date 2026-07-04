#!/usr/bin/env python3
# 运行方法：cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && conda run -n regular_bioinfo python 1_produce_nf-core_Samplesheet.py
import pandas as pd
import glob
import os

fastq_base_dir = "/Work_bio/dropbox/Dropbox_Data/Quote_06102601_Yue_Liu/"
original_sheet = "../Sample_Sheet_Yue_Liu.xlsx"
output_samplesheet = "nf_core_samplesheet.csv"

df = pd.read_excel(original_sheet)
df = df[df['Group'].notna()].reset_index(drop=True)
# Replicate number = running count within each Group (matches raw FASTQ naming convention)
df['rep'] = df.groupby('Group').cumcount() + 1

data = []
for _, row in df.iterrows():
    group = str(row['Group']).strip()
    rep = row['rep']
    sample_id = f"{group}_{rep}"

    pattern_r1 = os.path.join(fastq_base_dir, f"*_{group}_{rep}_S*_L*_R1_001.fastq.gz")
    pattern_r2 = os.path.join(fastq_base_dir, f"*_{group}_{rep}_S*_L*_R2_001.fastq.gz")
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

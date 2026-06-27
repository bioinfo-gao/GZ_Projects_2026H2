#!/usr/bin/env python3
# 运行方法：cd /home/gao/projects_2026H2/6_jinlong_mouse/scripts && python 1_produce_nf-core_Samplesheet.py
import pandas as pd
import glob
import os

fastq_base_dir = "/home/gao/Dropbox/Jinlong/"
original_csv   = "../jinlong.csv"
output_samplesheet = "nf_core_samplesheet.csv"

df = pd.read_csv(original_csv)
samples = df['Name in File'].unique()

def name_to_dir(name):
    """Map 'Name in File' (e.g. 902 / A) to actual directory name (J_902 / A)."""
    s = str(name).strip()
    if s.isdigit():
        return f"J_{s}"
    return s

data = []
for sample in samples:
    dir_name  = name_to_dir(sample)
    sample_dir = os.path.join(fastq_base_dir, dir_name)

    if not os.path.exists(sample_dir):
        print(f"WARNING: directory not found: {sample_dir}")
        continue

    r1_files = sorted(glob.glob(os.path.join(sample_dir, f"{dir_name}_*_1.fq.gz")))
    r2_files = sorted(glob.glob(os.path.join(sample_dir, f"{dir_name}_*_2.fq.gz")))

    if r1_files and r2_files:
        data.append([dir_name, r1_files[0], r2_files[0], 'reverse'])
    else:
        print(f"WARNING: FASTQ files not found for {dir_name} in {sample_dir}")

final_df = pd.DataFrame(data, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
final_df.to_csv(output_samplesheet, index=False)
print(f"Samplesheet saved: {output_samplesheet}")
print(f"Total samples: {len(data)}")
print(final_df.to_string())

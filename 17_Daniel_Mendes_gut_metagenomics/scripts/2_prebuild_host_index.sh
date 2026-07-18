#!/usr/bin/env bash
# 预建 mouse GRCm39 Bowtie2 去宿主索引（taxprofiler 短读长去宿主用），存共享参考目录长期复用，
# 避免每个项目让 pipeline 现场重建。建成后 taxprofiler 用 --shortread_hostremoval_index 指向该目录。
set -euo pipefail
REF=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa
IDXDIR=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/bowtie2_index
mkdir -p "$IDXDIR"
if ls "$IDXDIR"/GRCm39.*.bt2l >/dev/null 2>&1 || ls "$IDXDIR"/GRCm39.*.bt2 >/dev/null 2>&1; then
    echo "[$(date)] Bowtie2 index already present in $IDXDIR — skip"; exit 0
fi
echo "[$(date)] building Bowtie2 index -> $IDXDIR/GRCm39"
conda run -n regular_bioinfo bowtie2-build --threads 20 "$REF" "$IDXDIR/GRCm39"
echo "[$(date)] done: $(ls -la $IDXDIR | grep -c bt2) index files"

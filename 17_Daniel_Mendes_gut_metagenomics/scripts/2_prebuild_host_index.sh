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
# regular_bioinfo 无 bowtie2；mag_biobakery 有但是 2.2.3(太老，proj18 已踩坑)。
# 改用 taxprofiler 实跑时用过的 singularity 容器 bowtie2 2.4.2。
SIF=/Work_bio/gao/configs/.singularity/depot.galaxyproject.org-singularity-bowtie2-2.4.2--py38h1c8e9b9_1.img
singularity exec --bind /Work_bio:/Work_bio "$SIF" bowtie2-build --threads 20 "$REF" "$IDXDIR/GRCm39"
echo "[$(date)] done: $(ls -la $IDXDIR | grep -c bt2) index files"

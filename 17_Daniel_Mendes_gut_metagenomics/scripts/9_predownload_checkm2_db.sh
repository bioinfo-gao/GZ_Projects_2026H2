#!/usr/bin/env bash
# Phase 2 预备 — 下载 CheckM2 DIAMOND 库(uniref100.KO.*.dmnd)到共享参考目录，
# 让 nf-core/mag 的 --run_checkm2 能出 completeness/contamination(否则只有 BUSCO)。
# CheckM2 不在 mag_biobakery env，故直接从 Zenodo 取 tar、解压出 .dmnd 放到目标目录。
# 目标：/Work_bio/references/Metagenomics/checkm2/<something>.dmnd (脚本 4 检测 *.dmnd)
set -uo pipefail
DEST=/Work_bio/references/Metagenomics/checkm2
TMP=/Work_bio/references/Metagenomics/checkm2/_dl_tmp
URL=https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
mkdir -p "$DEST" "$TMP"

# 已就位则跳过（幂等）
if ls "$DEST"/*.dmnd >/dev/null 2>&1; then
  echo "[checkm2] .dmnd already present: $(ls "$DEST"/*.dmnd)"; exit 0
fi

echo "[$(date)] [checkm2] downloading DB tar ..."
wget -c -q --show-progress -O "$TMP/checkm2_database.tar.gz" "$URL" || { echo "[checkm2] DOWNLOAD FAIL"; exit 1; }
echo "[$(date)] [checkm2] extracting ..."
tar -xzf "$TMP/checkm2_database.tar.gz" -C "$TMP" || { echo "[checkm2] EXTRACT FAIL"; exit 1; }

# tar 内路径形如 CheckM2_database/uniref100.KO.1.dmnd —— 把 .dmnd 移到 DEST 顶层
DMND=$(find "$TMP" -name "*.dmnd" | head -1)
if [ -z "$DMND" ]; then echo "[checkm2] no .dmnd found in tar"; exit 1; fi
mv "$DMND" "$DEST"/
rm -rf "$TMP"
echo "[$(date)] [checkm2] DONE -> $(ls -lh "$DEST"/*.dmnd | awk '{print $5, $NF}')"
echo "CHECKM2_DB_READY"

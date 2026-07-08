# nf-core/mag 2.5.1 — 安装与测试记录（worked example for `/tax-resemb-mag`）

Date: 2026-07-08 · Prepared by: Zhen Gao, PhD, Athenomics · Platform: Linux HPC server

This folder is the validated worked example the `/tax-resemb-mag` skill references.
It records exactly how nf-core/mag 2.5.1 was installed into the `mag_biobakery`
conda env and validated with the built-in test profile.

---

## 1. What was installed

| Item | Detail |
| :--- | :---: |
| nextflow + openjdk | `mamba install -n mag_biobakery -c bioconda -c conda-forge nextflow` → **nextflow 26.04.4 + openjdk 23** |
| Pipeline | `nextflow pull nf-core/mag -r 2.5.1` (cached in `~/.nextflow/assets/nf-core/mag`, revision e72890089a) |
| Container engine | system apptainer/singularity 1.4.x + docker; run with `-profile singularity` |
| Shared image cache | `/Work_bio/gao/configs/.singularity` (shared across sarek/taxprofiler/mag) |
| BUSCO lineage (local) | `/Work_bio/references/Metagenomics/busco/bacteria_odb10.2024-01-08.tar.gz` |

**PATH caveat:** `mag_biobakery`'s own nextflow (26.04.4) is shadowed by
`regular_bioinfo/bin/nextflow` (25.10.4) in PATH. Scripts MUST call the absolute
path `/Work_bio/gao/configs/.conda/envs/mag_biobakery/bin/nextflow` and set
`JAVA_HOME=/Work_bio/gao/configs/.conda/envs/mag_biobakery`.

---

## 2. Three real problems hit during validation (all fixed)

1. **Config parsing failed on nextflow 26.04.4.** mag 2.5.1 uses the legacy nf-core
   template (`def check_max(obj, type)` in nextflow.config). Nextflow 26.04.x defaults
   to the new strict v2 config parser, which rejects it:
   `ERROR ~ Config parsing failed ... Unexpected input: '('`.
   **Fix:** `export NXF_SYNTAX_PARSER=v1` (classic parser). This is the single most
   important difference vs taxprofiler 2.0.1 (which uses the new template).

2. **Stale BUSCO test-DB URL.** The `test` profile hardcodes
   `bacteria_odb10.2020-03-06`, which busco-data.ezlab.org now 301-redirects to S3
   and 404s → workflow init dies with `No such file or directory`.
   **Fix:** downloaded the current `bacteria_odb10.2024-01-08.tar.gz` locally and
   passed `--busco_db <local tar.gz>`.

3. **CENTRIFUGE `/tmp` FIFO collision (exit 17).** Concurrent centrifuge containers
   each create a decompression FIFO named `/tmp/<in-container-PID>.inpipe1` on the
   shared host `/tmp`; in-container PIDs repeat across containers, so the second task
   dies with `mkfifo(/tmp/NN.inpipe1) failed`.
   **Fix:** `singularity.runOptions = '--writable-tmpfs'` (private writable /tmp per
   container). After this, centrifuge produced output for both test samples.

---

## 3. Validation result (`-profile test,singularity`)

```
-[nf-core/mag] Pipeline completed successfully, but with errored process(es)-
Duration : 7m 55s   Succeeded : 157   Ignored : 1   Failed : 1
```
- **Failed : 1** = `MAXBIN2 (test_minigut_sample2)` exit 255 — benign: MaxBin2 errors on
  the tiny synthetic contig set, and mag sets `errorStrategy ignore` for it. Not a real
  failure; pipeline exit code was 0.
- Validated end-to-end stages: fastp QC → Kraken2 + Centrifuge (both samples) →
  MEGAHIT + SPAdes assembly + QUAST → Prodigal → MetaBAT2 binning → BUSCO bin QC →
  BIN_SUMMARY → Prokka annotation → MultiQC.
- GTDB-Tk skipped by the test profile (`--skip_gtdbtk`, ~100GB db); CONCOCT skipped
  (`--skip_concoct`); Krona skipped.

Output tree under `test_run/results/`: `QC_shortreads/ Assembly/ GenomeBinning/
Taxonomy/{kraken2,centrifuge} Annotation/ multiqc/ pipeline_info/`.

---

## 4. Files in this folder

```
15_mag_setup/
├── scripts/
│   ├── local_resources.config   ← singularity(--writable-tmpfs) + 24c/96G caps
│   └── run_mag_test.sh          ← the validated test-run launcher (tmux + v1 parser + local busco_db)
├── docs/mag_setup_record_0708.md (this file)
└── test_run/                    ← test-profile run (results/ + work/ + logs)
```

The reusable operational knowledge is captured in the `/tax-resemb-mag` skill
(`~/.claude/commands/tax-resemb-mag.md`).

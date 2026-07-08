# nf-core/mag 5.4.2 — 安装与测试记录（worked example for `/tax-resemb-mag`）

Date: 2026-07-08 · Prepared by: Zhen Gao, PhD, Athenomics · Platform: Linux HPC server

This folder is the validated worked example the `/tax-resemb-mag` skill references.
It records how nf-core/mag **5.4.2** (the current release) was set up in the
`mag_biobakery` conda env and validated with the built-in test profile.

> This replaces an earlier throwaway attempt on mag 2.5.1 (legacy template), which was
> too old; 5.4.2 is the modern-template pipeline and needs none of 2.5.1's hacks.

---

## 1. What was installed

| Item | Detail |
| :--- | :---: |
| nextflow + openjdk | `mamba install -n mag_biobakery -c bioconda -c conda-forge nextflow` → **nextflow 26.04.4 + openjdk 23** |
| Pipeline | `nextflow pull nf-core/mag -r 5.4.2` (cached in `~/.nextflow/assets/nf-core/mag`, revision 5dabb0159a) |
| Min Nextflow required | mag 5.4.2 declares `nextflowVersion = '!>=25.04.2'` — 26.04.4 satisfies it |
| Container engine | system apptainer/singularity 1.4.x + docker; run with `-profile singularity` |
| Shared image cache | `/Work_bio/gao/configs/.singularity` (shared across sarek/taxprofiler/mag) |

**PATH caveat:** `mag_biobakery`'s own nextflow (26.04.4) is shadowed by
`regular_bioinfo/bin/nextflow` (25.10.4) in PATH. Scripts MUST call the absolute path
`/Work_bio/gao/configs/.conda/envs/mag_biobakery/bin/nextflow` and set
`JAVA_HOME=/Work_bio/gao/configs/.conda/envs/mag_biobakery`.

---

## 2. Why 5.4.2 needs none of the 2.5.1 workarounds

The earlier 2.5.1 setup required three hacks. All are obsolete in 5.4.2:

| 2.5.1 problem (legacy template) | 5.4.2 status |
| :--- | :---: |
| `NXF_SYNTAX_PARSER=v1` (legacy `check_max()` broke under Nextflow 26.x v2 parser) | **Not needed** — 5.4.2 uses the modern template; config parses cleanly |
| `--max_cpus/--max_memory` | **Removed** — resource ceilings come only from `process.resourceLimits` |
| CENTRIFUGE `/tmp` FIFO collision (`--writable-tmpfs`) | **Moot** — mag 5.x removed built-in read-level Kraken2/Centrifuge classification (use `/taxnom` for community profiling) |
| Stale BUSCO test URL override | **Not needed** — the test profile bundles its own mini DBs (BUSCO, GTDB-Tk mockup, CAT) |

---

## 3. Validation result (`-profile test,singularity`)

```
-[nf-core/mag] Pipeline completed successfully, but with errored process(es)-
Duration : 18m 15s   Succeeded : 195   Ignored : 1   Failed : 1
```
- **Failed : 1** = `MAXBIN2 (test_minigut_sample2)` exit 255 — benign: MaxBin2 errors on the
  tiny synthetic contig set and mag sets `errorStrategy ignore`. Pipeline exit code was 0.
- `WARN: Access to undefined parameter 'skip_tiara'` also appears — benign (Tiara euk
  classification off by default).
- End-to-end stages validated: fastp QC → Bowtie2 host removal → MEGAHIT + metaSPAdes
  assembly + QUAST → Prodigal → **6 binners** (MetaBAT2, MaxBin2, CONCOCT, COMEBin,
  MetaBinner, SemiBin2; CONCOCT/COMEBin/MetaBinner skipped by the default test profile) →
  **BUSCO** bin QC → **GTDB-Tk** (mockup db) → **CAT** bin classification → Prokka → MultiQC.
- Key outputs confirmed on disk: `GenomeBinning/QC/busco_summary.tsv`,
  `Taxonomy/GTDB-Tk/gtdbtk_summary.tsv`, `Taxonomy/CAT/bat_summary.tsv`,
  plus `Assembly/ GenomeBinning/ Annotation/ QC_shortreads/ multiqc/`.

**Real-project note:** the test uses a few-MB GTDB-Tk *mockup* db. Real bin classification
needs the full **GTDB r226 (~102GB)** — download once to a shared references dir and pass
`--gtdb_db`, or `--skip_gtdbtk` until it's staged.

---

## 4. Files in this folder

```
15_mag_setup/
├── scripts/
│   ├── local_resources.config   ← singularity cache + process.resourceLimits (24c/96G)
│   └── run_mag_test.sh          ← the validated test launcher (tmux + absolute nextflow path)
├── docs/mag_setup_record_0708.md (this file)
└── test_run/                    ← test-profile run (results/ + logs; work/ pruned after validation)
```

Reusable operational knowledge lives in the `/tax-resemb-mag` skill
(`~/.claude/commands/tax-resemb-mag.md`).

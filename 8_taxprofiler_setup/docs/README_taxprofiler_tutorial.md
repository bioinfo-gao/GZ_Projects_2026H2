# nf-core/taxprofiler — Local Installation & Test-Run Tutorial

Pipeline source: https://nf-co.re/taxprofiler/2.0.1
GitHub: https://github.com/nf-core/taxprofiler (pinned revision `2.0.1`)

nf-core/taxprofiler is a shotgun metagenomics taxonomic classification/profiling
pipeline. It runs short- and/or long-read QC, host removal, then fans reads out
to any combination of classifiers/profilers (Kraken2, Bracken, KrakenUniq,
Centrifuge, MetaPhlAn, Kaiju, DIAMOND, MALT, mOTUs, ganon, KMCP, sylph, MELON,
MetaCache) and standardises/aggregates results (MultiQC, Krona, taxpasta).

## 1. Environment

A dedicated mamba env `taxprofiler` was created (separate from the daily-driver
`regular_bioinfo` env) so the Nextflow/nf-core-tools version used for this
pipeline can be upgraded independently:

```bash
mamba create -n taxprofiler -y nextflow=25.10.* nf-core=4.0.* openjdk=17
```

Installed (verified 2026-07-01):

| Tool | Version | Path |
| :--- | :---: | :---: |
| nextflow | 25.10.4 | `/home/gao/.conda/envs/taxprofiler/bin/nextflow` |
| nf-core tools | 4.0.2 | same env |
| openjdk | 23 (already satisfies `>=17`) | same env |

Actual bioinformatics tools (Kraken2, MetaPhlAn, etc.) are **not** installed
into this conda env — nf-core pipelines run each process inside its own
Singularity/Apptainer container, pulled automatically by Nextflow. The env
only needs Nextflow + Java to orchestrate the workflow.

Container engine: this server already has `apptainer 1.4.5` (aliased as
`singularity`) and `docker 29.3.0` available system-wide; we use the
`singularity` execution profile (no root/daemon required, matches the
convention already used for nf-core/sarek in project `4_wgs_human_immu`).

## 2. Project layout

```
8_taxprofiler_setup/
├── docs/README_taxprofiler_tutorial.md   ← this file
├── configs/local_resources.config        ← CPU/RAM caps + singularity cache dir
├── scripts/1_run_taxprofiler_test_profile.sh  ← reproducible test-run launcher
├── logs/taxprofiler_test_run.log         ← full run log (tee'd from tmux)
└── test_run/
    ├── testdata/database_v2.1_taxprofiler2.0.1.csv  ← corrected test database sheet
    ├── work/                             ← Nextflow work dir (intermediate files)
    └── outdir/                           ← pipeline outputs (multiqc, taxpasta, etc.)
```

## 3. Singularity image cache

Reuses the shared cache already set up for other nf-core pipelines on this
machine, so images already pulled for e.g. sarek aren't re-downloaded:

```
NXF_SINGULARITY_CACHEDIR=/Work_bio/gao/configs/.singularity
```

This must be `export`-ed **inside** the tmux command string, not in an outer
shell script — if the tmux server is already resident, a new session inherits
the server's original environment, not variables exported by a wrapper script
right before `tmux new`. (Same pitfall hit previously with nf-core/sarek.)

## 4. Resource limits (`configs/local_resources.config`)

Server: AMD Threadripper 2990WX, 32 physical cores / 64 threads, 125 GB RAM.
Working-environment policy caps total concurrent usage at **28 physical
cores / 56 threads** to leave headroom for OS/SSH/interactive work. This
config caps the local Nextflow executor and per-process ceiling well inside
that limit (24 cores / 96 GB), on top of the pipeline's own `test` profile
process limits (4 cpus / 15 GB / 6h per process, from `conf/test.config`):

```groovy
executor {
    cpus   = 24
    memory = '96.GB'
}
process {
    resourceLimits = [cpus: 24, memory: 96.GB, time: 12.h]
}
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/Work_bio/gao/configs/.singularity'
}
```

## 5. Known issue: pinned-release test profile vs. moving test-datasets branch

The pipeline's built-in `test` profile (`conf/test.config`) points
`--databases` at `nf-core/test-datasets`' `taxprofiler` branch, which is a
**moving target** shared with the pipeline's `dev` branch (not frozen per
release tag). At the time this was run (2026-07-01), that CSV contained a
`centrifuger` database row — `centrifuger` is a classifier only recognised
by taxprofiler ≥2.1.x's parameter schema. Running the *pinned* `2.0.1`
release against the *current* test-datasets HEAD fails schema validation:

```
ERROR ~ Validation of pipeline parameters failed!
* --databases: Entry 6: Error for field 'tool' (centrifuger): Expected any of
  [bracken, centrifuge, diamond, ganon, kaiju, kmcp, kraken2, krakenuniq,
   malt, metaphlan, motus, sylph, melon, metacache] (Invalid tool name...)
```

**Fix applied**: downloaded the same CSV and removed the one `centrifuger`
row, saved as `test_run/testdata/database_v2.1_taxprofiler2.0.1.csv`, and
pass it explicitly via `--databases` (overriding the test profile default).
All other test-profile parameters (fastq samplesheet, host-removal reference,
per-tool `run_*` flags, etc.) are left as pipeline defaults.

If/when a real project database sheet is built, this same `--databases
<path/to/csv>` mechanism is how you point taxprofiler at production Kraken2 /
MetaPhlAn / etc. databases — see column format below.

## 6. How to (re)run the test profile

```bash
bash /home/gao/projects_2026H2/8_taxprofiler_setup/scripts/1_run_taxprofiler_test_profile.sh
```

This launches, inside tmux session `taxprofiler_test`:

```bash
nextflow run nf-core/taxprofiler \
  -r 2.0.1 \
  -profile test,singularity \
  -c configs/local_resources.config \
  --databases test_run/testdata/database_v2.1_taxprofiler2.0.1.csv \
  --outdir test_run/outdir \
  -work-dir test_run/work \
  -resume
```

Check progress:
```bash
tmux attach -t taxprofiler_test        # live view (Ctrl-B D to detach)
tail -f logs/taxprofiler_test_run.log  # or just tail the log
```

`-resume` is included so a re-run after any interruption reuses already
completed/cached work instead of restarting from scratch.

## 7. Running on real samples

For a production run, the two required inputs are:

1. **`--input samplesheet.csv`** — columns: `sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta`
   (see `assets/samplesheet.csv` in the pipeline repo for the exact template).
2. **`--databases databases.csv`** — columns: `tool,db_name,db_params,db_type,db_path`
   (`db_type` is `short`, `long`, or `short;long`; `db_path` is a local path or
   tarball/URL to the pre-built reference database for that classifier).

Typical invocation once real data/databases are staged:

```bash
nextflow run nf-core/taxprofiler \
  -r 2.0.1 \
  -profile singularity \
  -c configs/local_resources.config \
  --input samplesheet.csv \
  --databases databases.csv \
  --outdir results \
  -work-dir work \
  -resume
```

Toggle which classifiers actually run with the `--run_<tool>` boolean flags
(e.g. `--run_kraken2 --run_metaphlan`) — only databases matching an enabled
tool are used. QC/host-removal toggles (`--perform_shortread_qc`,
`--perform_shortread_hostremoval`, `--hostremoval_reference`, etc.) follow
the same pattern; full parameter list: `nf-core launch nf-core/taxprofiler`
or the schema at `nextflow_schema.json` / https://nf-co.re/taxprofiler/2.0.1/parameters.

### 7.1 Verified parameter names (avoid silently-ignored flags)

On 2026-07-02, checked a batch of plausible-looking parameters with
`nextflow run ... -preview` (parses/prints params without executing). Unknown
parameter names are **silently dropped** by the nf-schema plugin — no error,
no warning, not even listed in the params summary — so a typo or a
carried-over flag from another pipeline (e.g. old-template sarek) looks like
it worked but does nothing. Findings:

| Flag tried | Valid? | Notes |
| :--- | :---: | :---: |
| `--perform_shortread_preprocessing` | ❌ doesn't exist | Silently dropped. Correct flag: **`--perform_shortread_qc`** (add `--perform_shortread_complexityfilter` for low-complexity filtering) |
| `--remove_host` | ❌ doesn't exist | Silently dropped. Correct flag: **`--perform_shortread_hostremoval`** / `--perform_longread_hostremoval`, and it **requires** `--hostremoval_reference <genome.fasta>` (or a pre-built index) to actually do anything — the boolean alone is a no-op |
| `--run_kraken2` | ✅ valid | All `run_*` toggles default off for a real (non-`test`) run; must set explicitly, and `databases.csv` needs a `tool=kraken2` row |
| `--run_bracken` | ✅ valid | Schema description says it "automatically triggers the required Kraken2 prerequisite step" — turning this on alone already implies Kraken2; combining with `--run_kraken2` is redundant but harmless |
| `--run_metaphlan` | ✅ valid | Needs a `tool=metaphlan` row in `databases.csv` |
| `--max_cpus 28` | ❌ invalid | taxprofiler 2.0.1 uses the **new nf-core template**, which dropped the old `check_max()` + `--max_cpus`/`--max_memory` mechanism (only used by old-template pipelines like sarek 3.8.1). `conf/base.config` here only has fixed `process.cpus`/`withLabel` blocks with no reference to `params.max_cpus` at all — the flag is a complete no-op |
| `--max_memory '100.GB'` | ❌ invalid | Same as above, no effect |

**Conclusion**: resource caps must be set via Nextflow's native
`process.resourceLimits` block — i.e. this project's `configs/local_resources.config`
(section 4) — not via `--max_cpus`/`--max_memory` CLI flags. **Reproduction
tip**: before adding a custom flag to any new-template nf-core pipeline run,
verify it with `nextflow run ... -preview` and confirm it actually shows up
in the printed params summary, rather than assuming it works because it did
on a different (possibly old-template) pipeline.

## 8. Result status

**Test run completed successfully at 2026-07-01 23:45:12** (after recovering
from one dataflow stall — see below): `Succeeded: 103, Cached: 76, Failed: 0`.
All 14 classifiers produced output; `test_run/outdir/multiqc/multiqc_report.html`
is the summary report. Full per-tool outputs, standardised taxpasta tables,
and `pipeline_info/` (execution report/timeline/DAG/trace for all 3 attempts)
are under `test_run/outdir/`.

### 8.1 Incident: one dataflow stall during the run, and how it was resolved

Partway through, the run silently stopped making progress after the two
`METAPHLAN_METAPHLAN` tasks completed at 23:12:22 — no new log lines, no new
work-dir files, near-idle load average, main Nextflow JVM thread parked on
`futex_wait_queue`. Confirmed as a genuine dataflow deadlock (not just slow)
via Nextflow's own internal log, which printed this every 5 minutes starting
at 23:16:17:

```
!! executor local > No more task to compute -- The following nodes are still active:
```

This is Nextflow's own admission that downstream aggregator processes
(`TAXPASTA_MERGE`, `MULTIQC`, the `STANDARDISATION_PROFILES` collectors, etc.)
were stuck waiting on channels that never closed, and the scheduler had no
runnable task — a known edge case in Nextflow's dataflow scheduler under the
large, branchy DAG this test profile produces (~14 classifiers). Not an
environment or resource-config issue.

Fix: `tmux kill-session -t taxprofiler_test`, then re-run the same
`-resume`-based launch script. ~30 already-completed tasks were reused
(`cached: N ✔` in the log) instead of recomputed, and the stuck branches
(Centrifuge, sylph) resumed normally and the pipeline finished in another
~6 minutes of net compute.

**Reproduction tip**: if a future run goes >15 minutes with no new log output,
check `scripts/.nextflow.log` for repeated `No more task to compute -- ...
still active` lines before concluding it's stuck — that message is the
authoritative signal. If seen, `tmux kill-session` + re-run the same
`-resume` command; no completed work is lost.

#!/bin/bash
# Purpose: validate the local nf-core/taxprofiler 2.0.1 installation, container
# pulls, and executor configuration using the pipeline's official built-in
# test profile (tiny public test fastq + tiny test databases for every
# supported classifier). No real metagenomic data or local databases needed —
# everything is fetched automatically from nf-core/test-datasets.

set -euo pipefail

PROJECT_DIR="/home/gao/projects_2026H2/8_taxprofiler_setup"
cd "$PROJECT_DIR/scripts"

tmux kill-session -t taxprofiler_test 2>/dev/null || true

# NOTE: export must happen INSIDE the tmux command string, not in this outer
# script. If the tmux server is already resident, a new session inherits the
# environment the server itself started with, not variables exported by an
# outer script right before `tmux new`. (Same pitfall documented for the
# nf-core/sarek setup in project 4_wgs_human_immu.)
#
# NOTE on --databases override: the pipeline's built-in test profile points
# --databases at nf-core/test-datasets' `taxprofiler` branch, which is a
# moving target shared with the pipeline's dev branch. As of this csv's
# current HEAD it contains a `centrifuger` row, a tool only supported from
# taxprofiler >=2.1.x — pinning to release 2.0.1 fails schema validation
# against it. We use a local copy of the same CSV with that one row removed
# (see docs/README_taxprofiler_tutorial.md for details).

tmux new -d -s taxprofiler_test "
  export NXF_OPTS='-Xms512m -Xmx2g';
  export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity';
  /home/gao/.conda/envs/taxprofiler/bin/nextflow run nf-core/taxprofiler \
    -r 2.0.1 \
    -profile test,singularity \
    -c $PROJECT_DIR/configs/local_resources.config \
    --databases $PROJECT_DIR/test_run/testdata/database_v2.1_taxprofiler2.0.1.csv \
    --outdir $PROJECT_DIR/test_run/outdir \
    -work-dir $PROJECT_DIR/test_run/work \
    -resume \
    2>&1 | tee $PROJECT_DIR/logs/taxprofiler_test_run.log
"

echo "Started nf-core/taxprofiler test run in tmux session 'taxprofiler_test'."
echo "Live log: tmux attach -t taxprofiler_test   (or tail -f $PROJECT_DIR/logs/taxprofiler_test_run.log)"

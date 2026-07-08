#!/usr/bin/env bash
# nf-core/mag 5.4.2 — built-in test-profile validation run.
# Self-relaunch into tmux; on first failure, auto-retry with -resume.
set -uo pipefail

PROJ=/home/gao/projects_2026H2/15_mag_setup
# Use mag_biobakery's OWN nextflow (26.04.4). PATH ordering otherwise shadows it with
# regular_bioinfo's 25.10.4 -- and mag 5.4.2 requires Nextflow >= 25.04.2, so pin explicitly.
NEXTFLOW=/Work_bio/gao/configs/.conda/envs/mag_biobakery/bin/nextflow
export JAVA_HOME=/Work_bio/gao/configs/.conda/envs/mag_biobakery

if [ -z "${TMUX:-}" ]; then
    tmux new-session -d -s mag_test "bash '${PROJ}/scripts/run_mag_test.sh' 2>&1 | tee '${PROJ}/test_run/mag_test.log'"
    echo "launched in tmux session: mag_test"
    exit 0
fi

# Export INSIDE the tmux command context: a long-lived tmux server inherits its own
# startup env, so an outer export would not reach this session (verified pitfall).
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'

# NOTE: mag 5.4.2 is a MODERN-template pipeline, so NO NXF_SYNTAX_PARSER=v1 is needed
# (that hack was only for the legacy 2.x config). The test profile also bundles its own
# BUSCO / GTDB-Tk-mockup / CAT mini-DBs, so no --busco_db override is needed here either.

cd "${PROJ}/test_run"

run_mag() {
    "$NEXTFLOW" run nf-core/mag -r 5.4.2 \
        -profile test,singularity \
        -c "${PROJ}/scripts/local_resources.config" \
        --outdir "${PROJ}/test_run/results" \
        -work-dir "${PROJ}/test_run/work" "$@"
}

if run_mag;         then echo "MAG_TEST_SUCCESS"; exit 0; fi
echo "first attempt failed -- retrying with -resume"
if run_mag -resume; then echo "MAG_TEST_SUCCESS"; exit 0; fi
echo "MAG_TEST_FAILED"
exit 1

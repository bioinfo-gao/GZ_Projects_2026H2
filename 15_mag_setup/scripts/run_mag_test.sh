#!/usr/bin/env bash
# nf-core/mag 2.5.1 — built-in test-profile validation run.
# Self-relaunch into tmux; on first failure, auto-retry with -resume.
set -uo pipefail

PROJ=/home/gao/projects_2026H2/15_mag_setup
# Use mag_biobakery's OWN nextflow (PATH ordering otherwise shadows it with regular_bioinfo's).
NEXTFLOW=/Work_bio/gao/configs/.conda/envs/mag_biobakery/bin/nextflow
export JAVA_HOME=/Work_bio/gao/configs/.conda/envs/mag_biobakery

if [ -z "${TMUX:-}" ]; then
    tmux new-session -d -s mag_test "bash '${PROJ}/scripts/run_mag_test.sh' 2>&1 | tee '${PROJ}/test_run/mag_test.log'"
    echo "launched in tmux session: mag_test"
    exit 0
fi

# Must export INSIDE the tmux command context: a long-lived tmux server inherits its
# own startup env, so an outer export would not reach this session (verified taxprofiler pitfall).
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'

# nf-core/mag 2.5.1 is a LEGACY-template pipeline (uses def check_max() in nextflow.config).
# Nextflow 26.04.x defaults to the new strict (v2) config parser, which rejects that syntax
# ("Unexpected input: '('"). Force the classic v1 parser so the legacy config parses. (verified)
export NXF_SYNTAX_PARSER=v1

cd "${PROJ}/test_run"

# The test profile hardcodes a BUSCO db URL (bacteria_odb10.2020-03-06) that busco-data.ezlab.org
# no longer serves (301 -> S3 404). Override with a locally-downloaded current lineage tarball.
BUSCO_DB=/Work_bio/references/Metagenomics/busco/bacteria_odb10.2024-01-08.tar.gz

run_mag() {
    "$NEXTFLOW" run nf-core/mag -r 2.5.1 \
        -profile test,singularity \
        -c "${PROJ}/scripts/local_resources.config" \
        --busco_db "$BUSCO_DB" \
        --outdir "${PROJ}/test_run/results" \
        -work-dir "${PROJ}/test_run/work" "$@"
}

if run_mag;         then echo "MAG_TEST_SUCCESS"; exit 0; fi
echo "first attempt failed -- retrying with -resume"
if run_mag -resume; then echo "MAG_TEST_SUCCESS"; exit 0; fi
echo "MAG_TEST_FAILED"
exit 1

#!/usr/bin/env bash
# Project 16 — orchestrator: wait for sarek to finish, then auto-run the downstream chain
# (step 4 annotate -> step 5 rare/functional filter -> step 6 HLA -> step 8 origin check).
# Autonomous per user directive 2026-07-16: only a genuine blocking failure stops the chain.
# Launch in tmux: tmux new-session -d -s pan_down "bash scripts/9_run_downstream_when_ready.sh"
set -uo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
cd "$PROJ"
SLOG="$PROJ/logs/sarek_run.log"
NLOG="$PROJ/.nextflow.log"
DLOG="$PROJ/logs/downstream_orchestrator.log"
mkdir -p "$PROJ/logs"

say() { echo "[$(date +'%m-%d %H:%M:%S')] $*" | tee -a "$DLOG"; }

# ---------------- 1) wait for sarek to complete ----------------
say "orchestrator start; waiting for sarek to finish..."
while true; do
    if grep -q "=== sarek OK" "$SLOG" 2>/dev/null || \
       grep -qE "Pipeline completed successfully|Succeeded +: *[0-9]" "$NLOG" 2>/dev/null; then
        say "sarek COMPLETED (success marker found)."; break
    fi
    if grep -q "=== sarek FAILED" "$SLOG" 2>/dev/null || \
       grep -qiE "Pipeline completed with errors|Execution aborted" "$NLOG" 2>/dev/null; then
        say "!! sarek FAILED — STOPPING (blocking error, needs attention). Not running downstream."
        exit 1
    fi
    # liveness: if the run tmux is gone AND no success marker, treat as abnormal exit
    if ! tmux has-session -t pan_wgs 2>/dev/null; then
        sleep 20
        if ! grep -q "=== sarek OK" "$SLOG" 2>/dev/null; then
            say "!! pan_wgs session gone without success marker — STOPPING for inspection."
            exit 1
        fi
    fi
    sleep 120
done

# ---------------- 2) downstream chain (independent steps keep going on soft failure) ----------------
run_step() {
    local name="$1"; shift
    say ">>> START $name"
    if "$@" >> "$DLOG" 2>&1; then say ">>> OK    $name"; return 0
    else say ">>> FAIL  $name (exit $?) — continuing with independent steps"; return 1; fi
}

# step 4 -> 5 are dependent (5 needs 4's annotated VCFs)
if run_step "step4_annotate_gnomad_clinvar" bash scripts/4_annotate_gnomad_clinvar.sh; then
    run_step "step5_rare_functional_filter" conda run -n regular_bioinfo python scripts/5_rare_functional_filter.py
else
    say "skip step5 (depends on step4 output)."
fi
# step 6 (HLA) and step 8 (origin) are independent of 4/5
run_step "step6_hla_typing"  bash scripts/6_hla_typing.sh
run_step "step8_origin_check" bash scripts/8_origin_check.sh

say "downstream chain finished. Deliverables to assemble next: annotation_gnomad_clinvar/, prioritised_variants/, hla_typing/, origin_check/"
say "ALL_DONE"

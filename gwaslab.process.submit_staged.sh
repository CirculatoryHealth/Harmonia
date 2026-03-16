#!/usr/bin/env bash
#
# gwaslab.process.submit_staged.sh — submit a chained per-study SLURM pipeline
#
# For each active dataset in the config file this script submits one SLURM job
# per pipeline stage, chained with --dependency=afterok so that:
#   • if a stage fails, all subsequent stages for that study are automatically
#     cancelled by SLURM (DependencyNeverSatisfied)
#   • other studies are completely independent and continue running
#
# Stage chain (per study):
#   preprocess → process-normalize → process-check-ref → process-infer-strand
#   → [process-assign-rsid, if --dbsnp is set in WORKER_FLAGS]
#   → process-check-af → qc → cojo
#
# Two resource tiers are read from gwas_list.txt (COL8/COL9 and COL10/COL11):
#
#   HEAVY tier  (COL8=MEM, COL9=TIME)   — large VCF sweeps:
#     process-infer-strand   ← 1KG VCF sweep
#     process-assign-rsid    ← dbSNP VCF sweep (largest)
#     process-check-af       ← 1KG VCF sweep
#
#   LIGHT tier  (COL10=MEM_LIGHT, COL11=TIME_LIGHT)  — moderate steps:
#     process-normalize      ← DataFrame ops only
#     process-check-ref      ← FASTA random access; can spike for wide/complex files
#     qc                     ← filtering + plots
#
#   Fixed (trivial, not configurable per study):
#     preprocess             : 32G  / 12:00:00  ← CSV load + standardise
#     cojo                   : 16G  /  4:00:00  ← file write only
#
# Script-level fallback defaults (used when COL8–COL11 are absent from config):
#   LIGHT fallback : 64G  / 24:00:00
#   HEAVY fallback : 128G / 96:00:00
#
# Set MEM_LIGHT higher for studies with many columns or complex allele structure
# (e.g. multi-ancestry meta-analyses may need 128G even at process-check-ref).
#
# Usage:
#   bash gwaslab.process.submit_staged.sh gwas_list.txt
#   bash gwaslab.process.submit_staged.sh gwas_list.txt --partition=highmem
#
# Any extra arguments after the config file are forwarded to every sbatch call
# (e.g. --partition, --account, --reservation).
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION — adjust for your HPC environment
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/gwaslab.process.array_for_submit.sh"
LOG_BASE="/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo"

# ── Submission log (tee stdout+stderr to a timestamped file) ──────────────────
SUBMIT_LOG="${LOG_BASE}/gwaslab.process.submit_staged_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${LOG_BASE}"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "Submission log: ${SUBMIT_LOG}"

# Flags always passed to the Python worker (must match across all stages).
# --dbsnp controls whether process-assign-rsid is submitted.
WORKER_FLAGS="--liftover --figures --threads 8 --dbsnp --qc --cojo --cojo-pos --cojo-id rsid --leads --fill-eaf"

# ── Fixed (trivial) stage resources — not per-study configurable ──────────────
MEM_PREPROCESS="32G";  TIME_PREPROCESS="12:00:00"  # CSV load + standardise only
MEM_COJO="16G";        TIME_COJO="4:00:00"          # file write only

# ── Script-level fallback defaults (used when config COL10/COL11 are absent) ─
MEM_LIGHT_DEFAULT="64G";   TIME_LIGHT_DEFAULT="24:00:00"  # normalize, check-ref, qc
MEM_HEAVY_DEFAULT="128G";  TIME_HEAVY_DEFAULT="96:00:00"  # infer-strand, assign-rsid, check-af

# Whether process-assign-rsid is included in the chain.
# Auto-detected from WORKER_FLAGS; set to 1 when --dbsnp is present.
USE_DBSNP=0
if [[ "${WORKER_FLAGS}" == *"--dbsnp"* ]]; then USE_DBSNP=1; fi

# ─────────────────────────────────────────────────────────────────────────────
# Arguments
# ─────────────────────────────────────────────────────────────────────────────
CONFIG="${1:?Usage: bash gwaslab.process.submit_staged.sh <config.txt> [extra sbatch args]}"
shift   # remaining args forwarded to every sbatch call

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config file not found: ${CONFIG}" >&2; exit 1
fi
if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "ERROR: worker script not found: ${WORKER_SCRIPT}" >&2; exit 1
fi

# ── Read valid lines ──────────────────────────────────────────────────────────
mapfile -t LINES < <(grep -v '^\s*#' "${CONFIG}" | grep -v '^\s*$')
NTOTAL="${#LINES[@]}"
if [[ "${NTOTAL}" -eq 0 ]]; then
    echo "ERROR: no valid entries found in ${CONFIG}" >&2; exit 1
fi

echo "Config      : ${CONFIG}"
echo "Studies     : ${NTOTAL}"
echo "Worker      : ${WORKER_SCRIPT}"
echo "Use dbSNP   : ${USE_DBSNP}"
if [[ $# -gt 0 ]]; then echo "Extra sbatch: $*"; fi
echo "──────────────────────────────────────────────────────────────────────────────"
printf "%-22s  %-8s %-10s  %-8s %-10s  %s\n" \
    "GWAS" "MEM_L" "TIME_L" "MEM_H" "TIME_H" "JOB CHAIN"
echo "──────────────────────────────────────────────────────────────────────────────"

SUBMITTED=0
SKIPPED=0

for LINE in "${LINES[@]}"; do

    # ── Parse config line ─────────────────────────────────────────────────────
    # COL1–COL9 are required; COL10–COL11 (MEM_LIGHT, TIME_LIGHT) are optional
    # and fall back to script defaults when absent.
    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS \
        MEM TIME MEM_LIGHT TIME_LIGHT \
        <<< "${LINE}"

    if [[ -z "${INPUT_PATH}" || -z "${GWAS_NAME}" || -z "${POPULATION}" || \
          -z "${BUILD}"      || -z "${MEM}"        || -z "${TIME}" ]]; then
        echo "SKIP (malformed — missing required field in COL1–COL9): ${LINE}" >&2
        (( SKIPPED++ )) || true
        continue
    fi

    # ── Resolve per-study resources, falling back to script defaults ──────────
    # LIGHT tier: COL10/COL11 → fallback to MEM_LIGHT_DEFAULT/TIME_LIGHT_DEFAULT
    _MEM_LIGHT="${MEM_LIGHT:-${MEM_LIGHT_DEFAULT}}"
    _TIME_LIGHT="${TIME_LIGHT:-${TIME_LIGHT_DEFAULT}}"

    # HEAVY tier: COL8/COL9 → fallback to MEM_HEAVY_DEFAULT/TIME_HEAVY_DEFAULT
    _MEM_HEAVY="${MEM:-${MEM_HEAVY_DEFAULT}}"
    _TIME_HEAVY="${TIME:-${TIME_HEAVY_DEFAULT}}"

    # ── Submit the chain ──────────────────────────────────────────────────────
    # The worker script (array_for_submit.sh) handles conda activation, path
    # definitions, and building the python command.  We pass LINE as $1 and the
    # stage name as $2.  Each sbatch captures the job ID for the next dependency.

    # preprocess — fixed resources (trivial: CSV load + column standardisation)
    JID_PRE=$(sbatch \
        --job-name="gl_${GWAS_NAME}_preprocess" \
        --mem="${MEM_PREPROCESS}" --time="${TIME_PREPROCESS}" \
        --output="${LOG_BASE}/${GWAS_NAME}_1_preprocess_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_1_preprocess_%j.err" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "preprocess" \
        | awk '{print $NF}')

    # process-normalize — LIGHT tier
    JID_NRM=$(sbatch \
        --job-name="gl_${GWAS_NAME}_normalize" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --output="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.err" \
        --dependency="afterok:${JID_PRE}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-normalize" \
        | awk '{print $NF}')

    # process-check-ref — LIGHT tier (can spike for wide/multi-ancestry files)
    JID_CHR=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkref" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --output="${LOG_BASE}/${GWAS_NAME}_3_checkref_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_3_checkref_%j.err" \
        --dependency="afterok:${JID_NRM}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-check-ref" \
        | awk '{print $NF}')

    # process-infer-strand — HEAVY tier (1KG VCF full sweep)
    JID_IST=$(sbatch \
        --job-name="gl_${GWAS_NAME}_inferstrand" \
        --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
        --output="${LOG_BASE}/${GWAS_NAME}_4_inferstrand_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_4_inferstrand_%j.err" \
        --dependency="afterok:${JID_CHR}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-infer-strand" \
        | awk '{print $NF}')

    # process-assign-rsid — HEAVY tier (dbSNP VCF full sweep; only if --dbsnp)
    if [[ "${USE_DBSNP}" -eq 1 ]]; then
        JID_RSI=$(sbatch \
            --job-name="gl_${GWAS_NAME}_assignrsid" \
            --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
            --output="${LOG_BASE}/${GWAS_NAME}_5_assignrsid_%j.out" \
            --error="${LOG_BASE}/${GWAS_NAME}_5_assignrsid_%j.err" \
            --dependency="afterok:${JID_IST}" \
            "$@" \
            "${WORKER_SCRIPT}" "${LINE}" "process-assign-rsid" \
            | awk '{print $NF}')
        JID_PREV_CHECKAF="${JID_RSI}"
    else
        JID_RSI="(skipped)"
        JID_PREV_CHECKAF="${JID_IST}"
    fi

    # process-check-af — HEAVY tier (1KG VCF full sweep)
    JID_CAF=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkaf" \
        --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
        --output="${LOG_BASE}/${GWAS_NAME}_6_checkaf_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_6_checkaf_%j.err" \
        --dependency="afterok:${JID_PREV_CHECKAF}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-check-af" \
        | awk '{print $NF}')

    # qc — LIGHT tier (QC filtering + plots; no VCF sweeps)
    JID_QC=$(sbatch \
        --job-name="gl_${GWAS_NAME}_qc" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --output="${LOG_BASE}/${GWAS_NAME}_7_qc_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_7_qc_%j.err" \
        --dependency="afterok:${JID_CAF}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "qc" \
        | awk '{print $NF}')

    # cojo — fixed resources (trivial: file write only)
    JID_COJO=$(sbatch \
        --job-name="gl_${GWAS_NAME}_cojo" \
        --mem="${MEM_COJO}" --time="${TIME_COJO}" \
        --output="${LOG_BASE}/${GWAS_NAME}_8_cojo_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_8_cojo_%j.err" \
        --dependency="afterok:${JID_QC}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "cojo" \
        | awk '{print $NF}')

    # ── Report ────────────────────────────────────────────────────────────────
    printf "%-22s  %-8s %-10s  %-8s %-10s  %s → %s → %s → %s → %s → %s → %s → %s\n" \
        "${GWAS_NAME}" \
        "${_MEM_LIGHT}" "${_TIME_LIGHT}" \
        "${_MEM_HEAVY}" "${_TIME_HEAVY}" \
        "${JID_PRE}" "${JID_NRM}" "${JID_CHR}" "${JID_IST}" \
        "${JID_RSI}" "${JID_CAF}" "${JID_QC}" "${JID_COJO}"

    (( SUBMITTED++ )) || true

done

echo "──────────────────────────────────────────────────────────────────────────────"
echo "Done: ${SUBMITTED} study chain(s) submitted, ${SKIPPED} line(s) skipped."
echo ""
echo "Monitor  :  squeue -u \$USER"
echo "Cancel   :  scancel --name=gl_<GWAS_NAME>_preprocess  (cancels whole chain)"
echo "Details  :  sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,MaxRSS"

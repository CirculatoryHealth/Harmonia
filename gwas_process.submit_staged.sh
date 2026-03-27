#!/usr/bin/env bash
#
# gwas_process.submit_staged.sh — submit a chained per-study SLURM pipeline
#
# For each active dataset in the config file this script submits SLURM jobs
# chained with --dependency=afterok so that:
#   • if a stage fails, all subsequent stages for that study are automatically
#     cancelled by SLURM (DependencyNeverSatisfied)
#   • other studies are completely independent and continue running
#
# ── Per-chromosome (chr-split) stage chain (default) ─────────────────────────
#
#   preprocess → process-normalize → process-split
#   → [SLURM array 1-26] process-check-ref
#   → [SLURM array 1-26] process-infer-strand
#   → [SLURM array 1-26] process-assign-rsid  (only if --dbsnp in WORKER_FLAGS)
#   → [SLURM array 1-26] process-check-af
#   → merge
#
#   Array task IDs: 1-22 (autosomes)  23=X  24=Y  25=nonPAR  26=MT
#   Tasks for chromosomes absent from a study exit gracefully (exit 0),
#   satisfying afterok for the next array stage.
#   afterok on an array job ID waits for ALL tasks in the array to succeed.
#
# ── Resource tiers ────────────────────────────────────────────────────────────
#
#   HEAVY tier  (COL8=MEM, COL9=TIME)   — per-chr VCF sweep jobs:
#     process-infer-strand   ← 1KG VCF sweep (per chr ~1/22 of full dataset)
#     process-assign-rsid    ← dbSNP VCF sweep (per chr ~1/22)
#     process-check-af       ← 1KG VCF sweep (per chr ~1/22)
#
#   LIGHT tier  (COL10=MEM_LIGHT, COL11=TIME_LIGHT)  — moderate steps:
#     process-normalize      ← DataFrame ops only
#     process-check-ref      ← FASTA random access (per chr)
#     merge                  ← concat + QC + plots (light — no VCF sweeps)
#
#   Fixed (trivial, not configurable per study):
#     preprocess             : 32G  / 00:30:00  ← CSV load + standardise
#     process-split          : 16G  / 00:30:00  ← parquet split only
#
# Script-level fallback defaults (used when COL8–COL11 are absent from config):
#   LIGHT fallback : 64G  / 24:00:00
#   HEAVY fallback : 128G / 96:00:00
#
# Set MEM_HEAVY lower than the old whole-genome tier since each per-chr job
# only processes ~1/22 of variants.  The defaults above are conservative
# upper bounds; tune down once you have baseline RSS measurements per study.
#
# Usage:
#   bash gwas_process.submit_staged.sh gwas_list.txt
#   bash gwas_process.submit_staged.sh gwas_list.txt --partition=highmem
#
# Any extra arguments after the config file are forwarded to every sbatch call
# (e.g. --partition, --account, --reservation).
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration — loaded from gwas2cojo.conf (next to this script).
# Copy gwas2cojo.conf.example → gwas2cojo.conf and fill in your paths once.
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/gwas_process.array_for_submit.sh"
CONF="${SCRIPT_DIR}/gwas2cojo.conf"
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: ${CONF} not found." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL
LOG_BASE="${OUT_BASE}"   # submit_staged.sh uses LOG_BASE for SLURM output paths
# Export the absolute conf path so the SLURM worker (array_for_submit.sh) can
# find it even after SLURM copies the script to its own spool directory.
export GWAS2COJO_CONF="${CONF}"

# ── Submission log (tee stdout+stderr to a timestamped file) ──────────────────
SUBMIT_LOG="${LOG_BASE}/gwas_process.submit_staged_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${LOG_BASE}"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "Submission log: ${SUBMIT_LOG}"

# Flags always passed to the Python worker (must match across all stages).
# --dbsnp controls whether process-assign-rsid is submitted.
WORKER_FLAGS="--liftover --figures --threads 8 --dbsnp --qc --cojo --cojo-pos --cojo-id rsid --ldsc --leads --fill-eaf"

# ── SLURM job settings ────────────────────────────────────────────────────────
NODES=1          # nodes per job (all stages are single-node)
CPUS=8           # CPUs per job — must match --threads N in WORKER_FLAGS above
MAIL_TYPE="FAIL" # NONE | BEGIN | END | FAIL | ALL
# EMAIL is loaded from gwas2cojo.conf

# ── Fixed (trivial) stage resources — not per-study configurable ──────────────
MEM_PREPROCESS="32G";  TIME_PREPROCESS="00:30:00"  # CSV load + standardise only
MEM_SPLIT="16G";       TIME_SPLIT="00:30:00"        # parquet split only

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
CONFIG="${1:?Usage: bash gwas_process.submit_staged.sh <config.txt> [extra sbatch args]}"
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
    # COL12 (EXTRA_FLAGS) is optional; '.' means no extra flags.  It is passed
    # through to the worker script (array_for_submit.sh) via LINE, which appends
    # any flags verbatim to the Python command.  Example: --keep-multiallelic
    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS \
        MEM TIME MEM_LIGHT TIME_LIGHT EXTRA_FLAGS \
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
    # stage name as $2.  When $SLURM_ARRAY_TASK_ID is set inside the worker,
    # it automatically appends --chrom $SLURM_ARRAY_TASK_ID to the Python call.
    #
    # Array job dependency semantics:
    #   --dependency=afterok:ARRAY_JOB_ID  waits for ALL tasks in the array.
    #   Tasks that exit 0 (chromosome absent) count as satisfied.

    # ── 1. preprocess — fixed resources (CSV load + column standardisation) ───
    JID_PRE=$(sbatch \
        --job-name="gl_${GWAS_NAME}_preprocess" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${MEM_PREPROCESS}" --time="${TIME_PREPROCESS}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_1_preprocess_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_1_preprocess_%j.err" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "preprocess" \
        | awk '{print $NF}')

    # ── 2. process-normalize — LIGHT tier ─────────────────────────────────────
    JID_NRM=$(sbatch \
        --job-name="gl_${GWAS_NAME}_normalize" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.err" \
        --dependency="afterok:${JID_PRE}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-normalize" \
        | awk '{print $NF}')

    # ── 3. process-split — fixed (trivial: parquet split) ─────────────────────
    JID_SPL=$(sbatch \
        --job-name="gl_${GWAS_NAME}_split" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${MEM_SPLIT}" --time="${TIME_SPLIT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_3_split_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_3_split_%j.err" \
        --dependency="afterok:${JID_NRM}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-split" \
        | awk '{print $NF}')

    # ── 4. process-check-ref — LIGHT tier, per-chr array ──────────────────────
    JID_CHR=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkref" \
        --array="1-26" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_4_checkref_%A_%a.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_4_checkref_%A_%a.err" \
        --dependency="afterok:${JID_SPL}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-check-ref" \
        | awk '{print $NF}')

    # ── 5. process-infer-strand — HEAVY tier, per-chr array ───────────────────
    JID_IST=$(sbatch \
        --job-name="gl_${GWAS_NAME}_inferstrand" \
        --array="1-26" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_5_inferstrand_%A_%a.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_5_inferstrand_%A_%a.err" \
        --dependency="afterok:${JID_CHR}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-infer-strand" \
        | awk '{print $NF}')

    # ── 6. process-assign-rsid — HEAVY tier, per-chr array (only if --dbsnp) ──
    if [[ "${USE_DBSNP}" -eq 1 ]]; then
        JID_RSI=$(sbatch \
            --job-name="gl_${GWAS_NAME}_assignrsid" \
            --array="1-26" \
            --nodes="${NODES}" --cpus-per-task="${CPUS}" \
            --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
            --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
            --output="${LOG_BASE}/${GWAS_NAME}_6_assignrsid_%A_%a.out" \
            --error="${LOG_BASE}/${GWAS_NAME}_6_assignrsid_%A_%a.err" \
            --dependency="afterok:${JID_IST}" \
            "$@" \
            "${WORKER_SCRIPT}" "${LINE}" "process-assign-rsid" \
            | awk '{print $NF}')
        JID_PREV_CHECKAF="${JID_RSI}"
    else
        JID_RSI="(skipped)"
        JID_PREV_CHECKAF="${JID_IST}"
    fi

    # ── 7. process-check-af — HEAVY tier, per-chr array ───────────────────────
    JID_CAF=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkaf" \
        --array="1-26" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_HEAVY}" --time="${_TIME_HEAVY}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_7_checkaf_%A_%a.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_7_checkaf_%A_%a.err" \
        --dependency="afterok:${JID_PREV_CHECKAF}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-check-af" \
        | awk '{print $NF}')

    # ── 8. merge — LIGHT tier (concat + QC + plots + leads + COJO) ────────────
    # afterok on an array job ID waits for ALL 26 tasks to succeed.
    JID_MRG=$(sbatch \
        --job-name="gl_${GWAS_NAME}_merge" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_8_merge_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_8_merge_%j.err" \
        --dependency="afterok:${JID_CAF}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "merge" \
        | awk '{print $NF}')

    # ── Report ────────────────────────────────────────────────────────────────
    printf "%-22s  %-8s %-10s  %-8s %-10s  pre=%s nrm=%s spl=%s chr=%s ist=%s rsi=%s caf=%s mrg=%s\n" \
        "${GWAS_NAME}" \
        "${_MEM_LIGHT}" "${_TIME_LIGHT}" \
        "${_MEM_HEAVY}" "${_TIME_HEAVY}" \
        "${JID_PRE}" "${JID_NRM}" "${JID_SPL}" "${JID_CHR}" \
        "${JID_IST}" "${JID_RSI}" "${JID_CAF}" "${JID_MRG}"

    (( SUBMITTED++ )) || true

done

echo "──────────────────────────────────────────────────────────────────────────────"
echo "Done: ${SUBMITTED} study chain(s) submitted, ${SKIPPED} line(s) skipped."
echo ""
echo "Monitor  :  squeue -u \$USER"
echo "           squeue -u \$USER -r                       (show per-task array status)"
echo "Details  :  sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,MaxRSS"
echo ""
echo "Cancel   :  scancel --name=gl_<GWAS_NAME>_preprocess (full chain for one study)"
echo "           scancel <ARRAY_JOB_ID>                    (all tasks in one array job)"
echo "           scancel <ARRAY_JOB_ID>_<TASK_ID>          (one chromosome, e.g. 987654_3)"
echo "           scancel <ARRAY_JOB_ID>_[1-5]              (range of chromosomes)"
echo "           scancel <ARRAY_JOB_ID>_[1,3,22]           (comma list of chromosomes)"
echo "──────────────────────────────────────────────────────────────────────────────"

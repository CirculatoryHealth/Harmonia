#!/usr/bin/env bash
#
# gwas_process.submit_staged_b37.sh — submit a chained per-study SLURM pipeline
#                                     with GRCh37/hg19 output
#
# Identical to gwas_process.submit_staged.sh except:
#   • WORKER_SCRIPT → gwas_process.array_for_submit_b37.sh
#   • Output lands in ${OUT_BASE}/b37/${GWAS_NAME}
#   • --output-build 19 is passed globally (via the worker) instead of --liftover:
#       BUILD=38 inputs → reverse-lifted to GRCh37/hg19
#       BUILD=19/37 inputs → passed through unchanged
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
# Usage:
#   bash gwas_process.submit_staged_b37.sh gwas_list_b37.txt
#   bash gwas_process.submit_staged_b37.sh gwas_list_b37.txt --partition=highmem
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/gwas_process.array_for_submit_b37.sh"
CONF="${SCRIPT_DIR}/gwas2cojo.conf"
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: ${CONF} not found." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL
LOG_BASE="${OUT_BASE}/b37"
export GWAS2COJO_CONF="${CONF}"

# ── Submission log ─────────────────────────────────────────────────────────────
SUBMIT_LOG="${LOG_BASE}/gwas_process.submit_staged_b37_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${LOG_BASE}"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "Submission log: ${SUBMIT_LOG}"

# Flags that describe what the b37 worker passes to every Python invocation.
# Note: --output-build 19 and --dbsnp are hardcoded in the worker script; this
# variable is used only to drive USE_DBSNP detection below.
WORKER_FLAGS="--output-build 19 --figures --threads 8 --dbsnp --qc --cojo --cojo-pos --cojo-id rsid --ldsc --leads --fill-eaf"

# ── SLURM job settings ────────────────────────────────────────────────────────
NODES=1
CPUS=8
MAIL_TYPE="FAIL"

# ── Fixed (trivial) stage resources ──────────────────────────────────────────
MEM_PREPROCESS="${MEM_PREPROCESS:-32G}";  TIME_PREPROCESS="${TIME_PREPROCESS:-04:00:00}"
MEM_SPLIT="${MEM_SPLIT:-16G}";            TIME_SPLIT="${TIME_SPLIT:-01:00:00}"

# ── Script-level fallback defaults ────────────────────────────────────────────
MEM_LIGHT_DEFAULT="64G";   TIME_LIGHT_DEFAULT="24:00:00"
MEM_HEAVY_DEFAULT="128G";  TIME_HEAVY_DEFAULT="96:00:00"

USE_DBSNP=0
if [[ "${WORKER_FLAGS}" == *"--dbsnp"* ]]; then USE_DBSNP=1; fi

# ─────────────────────────────────────────────────────────────────────────────
# Arguments
# ─────────────────────────────────────────────────────────────────────────────
CONFIG="${1:?Usage: bash gwas_process.submit_staged_b37.sh <config.txt> [extra sbatch args]}"
shift

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config file not found: ${CONFIG}" >&2; exit 1
fi
if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "ERROR: worker script not found: ${WORKER_SCRIPT}" >&2; exit 1
fi

mapfile -t LINES < <(grep -v '^\s*#' "${CONFIG}" | grep -v '^\s*$')
NTOTAL="${#LINES[@]}"
if [[ "${NTOTAL}" -eq 0 ]]; then
    echo "ERROR: no valid entries found in ${CONFIG}" >&2; exit 1
fi

echo "Config      : ${CONFIG}"
echo "Studies     : ${NTOTAL}"
echo "Worker      : ${WORKER_SCRIPT}"
echo "Output base : ${LOG_BASE}"
echo "Use dbSNP   : ${USE_DBSNP}"
if [[ $# -gt 0 ]]; then echo "Extra sbatch: $*"; fi
echo "──────────────────────────────────────────────────────────────────────────────"
printf "%-22s  %-8s %-10s  %-8s %-10s  %s\n" \
    "GWAS" "MEM_L" "TIME_L" "MEM_H" "TIME_H" "JOB CHAIN"
echo "──────────────────────────────────────────────────────────────────────────────"

SUBMITTED=0
SKIPPED=0

for LINE in "${LINES[@]}"; do

    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS \
        MEM TIME MEM_LIGHT TIME_LIGHT EXTRA_FLAGS \
        <<< "${LINE}"

    if [[ -z "${INPUT_PATH}" || -z "${GWAS_NAME}" || -z "${POPULATION}" || \
          -z "${BUILD}"      || -z "${MEM}"        || -z "${TIME}" ]]; then
        echo "SKIP (malformed — missing required field in COL1–COL9): ${LINE}" >&2
        (( SKIPPED++ )) || true
        continue
    fi

    _MEM_LIGHT="${MEM_LIGHT:-${MEM_LIGHT_DEFAULT}}"
    _TIME_LIGHT="${TIME_LIGHT:-${TIME_LIGHT_DEFAULT}}"
    _MEM_HEAVY="${MEM:-${MEM_HEAVY_DEFAULT}}"
    _TIME_HEAVY="${TIME:-${TIME_HEAVY_DEFAULT}}"

    # ── 1. preprocess ─────────────────────────────────────────────────────────
    JID_PRE=$(sbatch \
        --job-name="b37_${GWAS_NAME}_preprocess" \
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
        --job-name="b37_${GWAS_NAME}_normalize" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_2_normalize_%j.err" \
        --dependency="afterok:${JID_PRE}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "process-normalize" \
        | awk '{print $NF}')

    # ── 3. process-split — fixed ───────────────────────────────────────────────
    JID_SPL=$(sbatch \
        --job-name="b37_${GWAS_NAME}_split" \
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
        --job-name="b37_${GWAS_NAME}_checkref" \
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
        --job-name="b37_${GWAS_NAME}_inferstrand" \
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
            --job-name="b37_${GWAS_NAME}_assignrsid" \
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
        --job-name="b37_${GWAS_NAME}_checkaf" \
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

    # ── 8. merge — LIGHT tier ─────────────────────────────────────────────────
    JID_MRG=$(sbatch \
        --job-name="b37_${GWAS_NAME}_merge" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${_MEM_LIGHT}" --time="${_TIME_LIGHT}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_8_merge_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_8_merge_%j.err" \
        --dependency="afterok:${JID_CAF}" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" "merge" \
        | awk '{print $NF}')

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
echo "Cancel   :  scancel --name=b37_<GWAS_NAME>_preprocess (full chain for one study)"
echo "           scancel <ARRAY_JOB_ID>                    (all tasks in one array job)"
echo "──────────────────────────────────────────────────────────────────────────────"

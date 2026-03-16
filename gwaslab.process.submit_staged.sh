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
# Resource defaults (override in the USER CONFIGURATION block below):
#   preprocess          : 32G  / 12:00:00
#   process-normalize   : 64G  / 24:00:00
#   process-check-ref   : 64G  / 24:00:00
#   process-infer-strand: 128G / 48:00:00   ← 1KG VCF sweep
#   process-assign-rsid : 256G / 96:00:00   ← dbSNP VCF sweep (largest)
#   process-check-af    : 128G / 48:00:00   ← 1KG VCF sweep
#   qc                  :  64G / 24:00:00
#   cojo                :  16G /  4:00:00
#
# The MEM and TIME columns from gwas_list.txt are applied to the two heaviest
# stages (process-infer-strand and process-assign-rsid) as a per-study override
# of the defaults above.
#
# Usage:
#   bash gwaslab.process.submit_staged.sh gwas_list.txt
#   bash gwaslab.process.submit_staged.sh gwas_list.txt --partition=highmem
#
# Any extra arguments are forwarded to every sbatch call (e.g. --partition,
# --account, --reservation).
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION — adjust for your HPC environment
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/gwaslab.process.array_for_submit.sh"
LOG_BASE="/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo"

# Flags always passed to the Python worker (must match across all stages).
# --dbsnp controls whether process-assign-rsid is submitted.
WORKER_FLAGS="--liftover --figures --threads 8 --dbsnp --qc --cojo --cojo-pos --cojo-id rsid --leads --fill-eaf"

# ── Per-stage resource defaults ───────────────────────────────────────────────
MEM_PREPROCESS="32G";           TIME_PREPROCESS="12:00:00"
MEM_NORMALIZE="64G";            TIME_NORMALIZE="24:00:00"
MEM_CHECKREF="64G";             TIME_CHECKREF="24:00:00"
# infer-strand and assign-rsid use MEM/TIME from the config file (per study).
# The values below are fallback defaults used when the config fields are missing.
MEM_INFERSTRAND_DEFAULT="128G"; TIME_INFERSTRAND_DEFAULT="48:00:00"
MEM_ASSIGNRSID_DEFAULT="256G";  TIME_ASSIGNRSID_DEFAULT="96:00:00"
MEM_CHECKAF="128G";             TIME_CHECKAF="48:00:00"
MEM_QC="64G";                   TIME_QC="24:00:00"
MEM_COJO="16G";                 TIME_COJO="4:00:00"

# Whether process-assign-rsid is included in the chain.
# Set to 1 when --dbsnp appears in WORKER_FLAGS (auto-detected below).
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
echo "──────────────────────────────────────────────────────────────────────"
printf "%-22s  %-10s  %-10s  %s\n" "GWAS" "MEM(heavy)" "TIME(heavy)" "JOB CHAIN"
echo "──────────────────────────────────────────────────────────────────────"

SUBMITTED=0
SKIPPED=0

for LINE in "${LINES[@]}"; do

    # ── Parse config line ─────────────────────────────────────────────────────
    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS MEM TIME \
        <<< "${LINE}"

    if [[ -z "${INPUT_PATH}" || -z "${GWAS_NAME}" || -z "${POPULATION}" || \
          -z "${BUILD}"      || -z "${MEM}"        || -z "${TIME}" ]]; then
        echo "SKIP (malformed — missing required field): ${LINE}" >&2
        (( SKIPPED++ )) || true
        continue
    fi

    # ── Per-study heavy-stage resources (from config; fallback to defaults) ───
    MEM_INFERSTRAND="${MEM:-${MEM_INFERSTRAND_DEFAULT}}"
    TIME_INFERSTRAND="${TIME:-${TIME_INFERSTRAND_DEFAULT}}"
    MEM_ASSIGNRSID="${MEM:-${MEM_ASSIGNRSID_DEFAULT}}"
    TIME_ASSIGNRSID="${TIME:-${TIME_ASSIGNRSID_DEFAULT}}"

    # ── Build the common python command prefix (excluding --stage) ────────────
    GWAS_DIR=$(dirname  "${INPUT_PATH}")
    INPUT_FILE=$(basename "${INPUT_PATH}")

    # Collect optional N arguments
    N_ARGS=()
    N_SET=0; NCASES_SET=0; NCONTROLS_SET=0
    if [[ -n "${N:-}"          && "${N}"          != "." ]]; then N_ARGS+=(--n         "${N}");        N_SET=1;        fi
    if [[ -n "${N_CASES:-}"    && "${N_CASES}"    != "." ]]; then N_ARGS+=(--n-cases   "${N_CASES}");  NCASES_SET=1;   fi
    if [[ -n "${N_CONTROLS:-}" && "${N_CONTROLS}" != "." ]]; then N_ARGS+=(--n-controls "${N_CONTROLS}"); NCONTROLS_SET=1; fi
    if [[ "${N_SET}" -eq 1 && "${NCASES_SET}" -eq 1 && "${NCONTROLS_SET}" -eq 1 ]]; then
        N_ARGS+=(--force-n)
    fi

    # Base python invocation — --stage is appended per-job below
    BASE_CMD=(
        python "${PYTHON_SCRIPT:-/hpc/local/Rocky8/dhl_ec/software/gwas2cojo/gwaslab.process.py}"
        --gwas        "${GWAS_NAME}"
        --input       "${INPUT_FILE}"
        --directory   "${GWAS_DIR}"
        --ref         "${REF_DIR:-/hpc/dhl_ec/data/references/gwaslab}"
        --output      "${OUT_BASE}/${GWAS_NAME}"
        --population  "${POPULATION}"
        --build       "${BUILD}"
        ${WORKER_FLAGS}
        "${N_ARGS[@]}"
    )

    # ── Helper: submit one stage, return job ID ───────────────────────────────
    submit_stage() {
        local stage_name="$1" mem="$2" time="$3" dep="$4"
        local dep_flag=""
        [[ -n "${dep}" ]] && dep_flag="--dependency=afterok:${dep}"

        sbatch \
            --job-name="gl_${GWAS_NAME}_${stage_name}" \
            --mem="${mem}" \
            --time="${time}" \
            --output="${LOG_BASE}/${GWAS_NAME}_${stage_name}_%j.out" \
            --error="${LOG_BASE}/${GWAS_NAME}_${stage_name}_%j.err" \
            ${dep_flag} \
            "$@_EXTRA" \
            --wrap="${BASE_CMD[*]} --stage ${stage_name}" \
            | awk '{print $NF}'
    }

    # ── Submit the chain ──────────────────────────────────────────────────────
    # Each sbatch call captures the job ID and passes it as dependency to the next.

    JID_PRE=$(sbatch \
        --job-name="gl_${GWAS_NAME}_preprocess" \
        --mem="${MEM_PREPROCESS}" --time="${TIME_PREPROCESS}" \
        --output="${LOG_BASE}/${GWAS_NAME}_preprocess_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_preprocess_%j.err" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage preprocess" \
        | awk '{print $NF}')

    JID_NRM=$(sbatch \
        --job-name="gl_${GWAS_NAME}_normalize" \
        --mem="${MEM_NORMALIZE}" --time="${TIME_NORMALIZE}" \
        --output="${LOG_BASE}/${GWAS_NAME}_normalize_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_normalize_%j.err" \
        --dependency="afterok:${JID_PRE}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage process-normalize" \
        | awk '{print $NF}')

    JID_CHR=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkref" \
        --mem="${MEM_CHECKREF}" --time="${TIME_CHECKREF}" \
        --output="${LOG_BASE}/${GWAS_NAME}_checkref_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_checkref_%j.err" \
        --dependency="afterok:${JID_NRM}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage process-check-ref" \
        | awk '{print $NF}')

    JID_IST=$(sbatch \
        --job-name="gl_${GWAS_NAME}_inferstrand" \
        --mem="${MEM_INFERSTRAND}" --time="${TIME_INFERSTRAND}" \
        --output="${LOG_BASE}/${GWAS_NAME}_inferstrand_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_inferstrand_%j.err" \
        --dependency="afterok:${JID_CHR}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage process-infer-strand" \
        | awk '{print $NF}')

    # process-assign-rsid: only submitted when --dbsnp is in WORKER_FLAGS
    if [[ "${USE_DBSNP}" -eq 1 ]]; then
        JID_RSI=$(sbatch \
            --job-name="gl_${GWAS_NAME}_assignrsid" \
            --mem="${MEM_ASSIGNRSID}" --time="${TIME_ASSIGNRSID}" \
            --output="${LOG_BASE}/${GWAS_NAME}_assignrsid_%j.out" \
            --error="${LOG_BASE}/${GWAS_NAME}_assignrsid_%j.err" \
            --dependency="afterok:${JID_IST}" \
            "$@" \
            --wrap="${BASE_CMD[*]} --stage process-assign-rsid" \
            | awk '{print $NF}')
        JID_PREV_CHECKAF="${JID_RSI}"
    else
        JID_RSI="(skipped)"
        JID_PREV_CHECKAF="${JID_IST}"
    fi

    JID_CAF=$(sbatch \
        --job-name="gl_${GWAS_NAME}_checkaf" \
        --mem="${MEM_CHECKAF}" --time="${TIME_CHECKAF}" \
        --output="${LOG_BASE}/${GWAS_NAME}_checkaf_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_checkaf_%j.err" \
        --dependency="afterok:${JID_PREV_CHECKAF}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage process-check-af" \
        | awk '{print $NF}')

    JID_QC=$(sbatch \
        --job-name="gl_${GWAS_NAME}_qc" \
        --mem="${MEM_QC}" --time="${TIME_QC}" \
        --output="${LOG_BASE}/${GWAS_NAME}_qc_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_qc_%j.err" \
        --dependency="afterok:${JID_CAF}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage qc" \
        | awk '{print $NF}')

    JID_COJO=$(sbatch \
        --job-name="gl_${GWAS_NAME}_cojo" \
        --mem="${MEM_COJO}" --time="${TIME_COJO}" \
        --output="${LOG_BASE}/${GWAS_NAME}_cojo_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_cojo_%j.err" \
        --dependency="afterok:${JID_QC}" \
        "$@" \
        --wrap="${BASE_CMD[*]} --stage cojo" \
        | awk '{print $NF}')

    # ── Report ────────────────────────────────────────────────────────────────
    printf "%-22s  %-10s  %-10s  %s → %s → %s → %s → %s → %s → %s → %s\n" \
        "${GWAS_NAME}" "${MEM_INFERSTRAND}" "${TIME_INFERSTRAND}" \
        "${JID_PRE}" "${JID_NRM}" "${JID_CHR}" "${JID_IST}" \
        "${JID_RSI}" "${JID_CAF}" "${JID_QC}" "${JID_COJO}"

    (( SUBMITTED++ )) || true

done

echo "──────────────────────────────────────────────────────────────────────"
echo "Done: ${SUBMITTED} study chain(s) submitted, ${SKIPPED} line(s) skipped."
echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Cancel study:  scancel --name=gl_<GWAS_NAME>_preprocess  (cancels chain)"
echo "Job details:   sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,MaxRSS"

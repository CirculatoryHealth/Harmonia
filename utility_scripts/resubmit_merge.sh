#!/usr/bin/env bash
#
# resubmit_merge.sh — resubmit a single pipeline stage for studies in a gwas_list file
#
# Designed to be called from ANY working directory (not only the gwas2cojo source
# directory).  The script locates harmonia.conf via three methods (in order):
#   1. --conf <path>  argument
#   2. HARMONIA_CONF environment variable
#   3. harmonia.conf next to this script's parent directory (auto-detect)
#
# Usage
# ─────
#   bash resubmit_merge.sh [OPTIONS] [-- extra sbatch args]
#
# Options
# ───────
#   --conf       PATH   Path to harmonia.conf
#   --gwas-list  FILE   Input list (gwas_list.txt format).
#                       Default: gwas_list_resubmit_merge.txt (next to gwas_list.txt)
#   --stage      STAGE  Pipeline stage to resubmit (default: merge).
#                       Any valid --stage value for harmonia.py:
#                         merge | qc | cojo | preprocess | process-normalize | ...
#   --mem        MEM    SLURM memory override (default: from gwas_list COL10 → 64G)
#   --time       TIME   SLURM time override   (default: from gwas_list COL11 → 06:00:00)
#   --dry-run           Print sbatch commands without submitting
#
# Examples
# ────────
#   # Run merge for all studies in gwas_list_resubmit_merge.txt:
#   bash /path/to/resubmit_merge.sh
#
#   # Override conf location (when not in the gwas2cojo directory):
#   bash /path/to/resubmit_merge.sh --conf /hpc/local/gwas2cojo/harmonia.conf
#
#   # Full rerun from preprocess using a different list:
#   bash /path/to/resubmit_merge.sh \
#       --gwas-list gwas_list_resubmit_full.txt \
#       --stage preprocess
#
#   # Dry-run to verify job parameters before submitting:
#   bash /path/to/resubmit_merge.sh --dry-run
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Script location (for auto-detecting sibling files) ────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# gwas2cojo root is one level up from utility_scripts/
ROOTDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Argument defaults ─────────────────────────────────────────────────────────
CONF_ARG=""
GWAS_LIST_ARG=""
STAGE="merge"
MEM_OVERRIDE=""
TIME_OVERRIDE=""
DRY_RUN=0

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --conf)       CONF_ARG="$2";       shift 2 ;;
        --gwas-list)  GWAS_LIST_ARG="$2";  shift 2 ;;
        --stage)      STAGE="$2";          shift 2 ;;
        --mem)        MEM_OVERRIDE="$2";   shift 2 ;;
        --time)       TIME_OVERRIDE="$2";  shift 2 ;;
        --dry-run)    DRY_RUN=1;           shift   ;;
        --)           shift;               break   ;;
        -*) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
        *)  break ;;
    esac
done
# Remaining args are passed verbatim to every sbatch call (e.g. --partition=highmem)
EXTRA_SBATCH=("$@")

# ── Locate harmonia.conf ──────────────────────────────────────────────────────
if [[ -n "${CONF_ARG}" ]]; then
    CONF="${CONF_ARG}"
elif [[ -n "${HARMONIA_CONF:-}" && -f "${HARMONIA_CONF}" ]]; then
    CONF="${HARMONIA_CONF}"
elif [[ -f "${ROOTDIR}/harmonia.conf" ]]; then
    CONF="${ROOTDIR}/harmonia.conf"
else
    echo "ERROR: harmonia.conf not found." >&2
    echo "       Use --conf /path/to/harmonia.conf, or export HARMONIA_CONF=..." >&2
    exit 1
fi
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: harmonia.conf not found at: ${CONF}" >&2
    exit 1
fi

# shellcheck source=../harmonia.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL
export HARMONIA_CONF="${CONF}"   # propagate to SLURM worker

# Derive the worker script from PYTHON_SCRIPT (handles installs in any location)
HARMONIA_ROOTDIR="$(dirname "${PYTHON_SCRIPT}")"
WORKER_SCRIPT="${HARMONIA_ROOTDIR}/harmonia.array_for_submit.sh"

if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "ERROR: worker script not found: ${WORKER_SCRIPT}" >&2
    echo "       Check PYTHON_SCRIPT in ${CONF}" >&2
    exit 1
fi

# ── Locate gwas_list input file ───────────────────────────────────────────────
if [[ -n "${GWAS_LIST_ARG}" ]]; then
    GWAS_LIST="${GWAS_LIST_ARG}"
else
    # Default: gwas_list_resubmit_merge.txt next to gwas_list.txt in ROOTDIR
    GWAS_LIST="${HARMONIA_ROOTDIR}/gwas_list_resubmit_merge.txt"
fi

if [[ ! -f "${GWAS_LIST}" ]]; then
    echo "ERROR: gwas_list not found: ${GWAS_LIST}" >&2
    echo "       Use --gwas-list /path/to/your_list.txt" >&2
    exit 1
fi

# ── Resource defaults (per-stage) ─────────────────────────────────────────────
# For merge / qc / cojo: light resources are sufficient.
# Caller can override with --mem / --time.
MEM_DEFAULT="64G"
TIME_DEFAULT="06:00:00"

# ── Submission log ─────────────────────────────────────────────────────────────
LOG_BASE="${OUT_BASE}"
mkdir -p "${LOG_BASE}"
SUBMIT_LOG="${LOG_BASE}/resubmit_${STAGE}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1

echo "========================================================"
echo "resubmit_merge.sh"
echo "Conf         : ${CONF}"
echo "GWAS list    : ${GWAS_LIST}"
echo "Stage        : ${STAGE}"
echo "Worker       : ${WORKER_SCRIPT}"
echo "Log base     : ${LOG_BASE}"
if [[ "${DRY_RUN}" -eq 1 ]]; then echo "Mode         : DRY-RUN (no jobs submitted)"; fi
echo "========================================================"

# ── Read valid lines ───────────────────────────────────────────────────────────
mapfile -t LINES < <(grep -v '^\s*#' "${GWAS_LIST}" | grep -v '^\s*$')
NTOTAL="${#LINES[@]}"
if [[ "${NTOTAL}" -eq 0 ]]; then
    echo "ERROR: no valid (uncommented) entries found in ${GWAS_LIST}" >&2
    exit 1
fi
echo "Studies      : ${NTOTAL}"
echo "──────────────────────────────────────────────────────────────────────────────"

SUBMITTED=0
SKIPPED=0

for LINE in "${LINES[@]}"; do

    # ── Parse config columns ───────────────────────────────────────────────────
    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS \
        MEM TIME MEM_LIGHT TIME_LIGHT EXTRA_FLAGS \
        <<< "${LINE}"

    if [[ -z "${GWAS_NAME:-}" ]]; then
        echo "SKIP (malformed — GWAS_NAME empty): ${LINE}" >&2
        (( SKIPPED++ )) || true
        continue
    fi

    # ── Resolve memory and time ────────────────────────────────────────────────
    # Priority: --mem/--time CLI override → COL10/COL11 (MEM_LIGHT/TIME_LIGHT) → defaults
    if [[ -n "${MEM_OVERRIDE}" ]]; then
        _MEM="${MEM_OVERRIDE}"
    elif [[ -n "${MEM_LIGHT:-}" && "${MEM_LIGHT}" != "." ]]; then
        _MEM="${MEM_LIGHT}"
    else
        _MEM="${MEM_DEFAULT}"
    fi

    if [[ -n "${TIME_OVERRIDE}" ]]; then
        _TIME="${TIME_OVERRIDE}"
    elif [[ -n "${TIME_LIGHT:-}" && "${TIME_LIGHT}" != "." ]]; then
        _TIME="${TIME_LIGHT}"
    else
        _TIME="${TIME_DEFAULT}"
    fi

    # ── Build sbatch command ───────────────────────────────────────────────────
    SBATCH_CMD=(
        sbatch
        --job-name="gl_${GWAS_NAME}_${STAGE//process-/}_resubmit"
        --nodes=1 --cpus-per-task=8
        --mem="${_MEM}" --time="${_TIME}"
        --mail-type=FAIL --mail-user="${EMAIL}"
        --output="${LOG_BASE}/${GWAS_NAME}_resubmit_${STAGE//process-/}_%j.out"
        --error="${LOG_BASE}/${GWAS_NAME}_resubmit_${STAGE//process-/}_%j.err"
        "${EXTRA_SBATCH[@]+"${EXTRA_SBATCH[@]}"}"
        "${WORKER_SCRIPT}" "${LINE}" "${STAGE}"
    )

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "DRY-RUN  ${GWAS_NAME}  mem=${_MEM}  time=${_TIME}  stage=${STAGE}"
        echo "         ${SBATCH_CMD[*]}"
        (( SUBMITTED++ )) || true
        continue
    fi

    JID=$("${SBATCH_CMD[@]}" | awk '{print $NF}')
    echo "Submitted  ${GWAS_NAME}  (${STAGE})  mem=${_MEM}  time=${_TIME}  jobid=${JID}"
    (( SUBMITTED++ )) || true

done

echo "──────────────────────────────────────────────────────────────────────────────"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "DRY-RUN: ${SUBMITTED} job(s) would be submitted, ${SKIPPED} skipped."
else
    echo "Done: ${SUBMITTED} job(s) submitted, ${SKIPPED} skipped."
    echo ""
    echo "Monitor : squeue -u \$USER"
    echo "Logs    : ${LOG_BASE}/<GWAS_NAME>_resubmit_${STAGE//process-/}_<JOB_ID>.out"
fi
echo "Submission log: ${SUBMIT_LOG}"

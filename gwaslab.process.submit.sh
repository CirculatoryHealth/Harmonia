#!/usr/bin/env bash
#
# gwaslab.process.submit.sh — submit one independent SLURM job per GWAS dataset.
#
# Memory and time limits are read from the config file (COL8 and COL9), so each
# dataset gets exactly the resources it needs.  Jobs are never killed because a
# sibling in an array timed out.
#
# Usage:
#   bash gwaslab.process.submit.sh gwas_list.txt
#   bash gwaslab.process.submit.sh gwas_list.txt --partition=highmem
#
# Any extra arguments are forwarded to every sbatch call (e.g. --partition,
# --account).  Per-job --mem, --time, --job-name, --output and --error are
# always set from the config file and cannot be overridden this way.
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/gwaslab.process.array_for_submit.sh"

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration — loaded from gwas2cojo.conf (next to this script).
# Copy gwas2cojo.conf.example → gwas2cojo.conf and fill in your paths once.
# ─────────────────────────────────────────────────────────────────────────────
CONF="${SCRIPT_DIR}/gwas2cojo.conf"
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: ${CONF} not found." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL
LOG_BASE="${OUT_BASE}"   # submit.sh uses LOG_BASE for SLURM output paths
# Export the absolute conf path so the SLURM worker (array_for_submit.sh) can
# find it even after SLURM copies the script to its own spool directory.
export GWAS2COJO_CONF="${CONF}"

# ── SLURM job settings ────────────────────────────────────────────────────────
NODES=1
CPUS=8
MAIL_TYPE="FAIL"   # NONE | BEGIN | END | FAIL | ALL

# ── Arguments ─────────────────────────────────────────────────────────────────
CONFIG="${1:?Usage: bash gwaslab.process.submit.sh <config.txt> [extra sbatch args]}"
shift   # remaining args forwarded to every sbatch call

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config file not found: ${CONFIG}" >&2
    exit 1
fi

if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "ERROR: worker script not found: ${WORKER_SCRIPT}" >&2
    exit 1
fi

# ── Read valid lines (non-blank, non-comment) ─────────────────────────────────
mapfile -t LINES < <(grep -v '^\s*#' "${CONFIG}" | grep -v '^\s*$')

NTOTAL="${#LINES[@]}"
if [[ "${NTOTAL}" -eq 0 ]]; then
    echo "ERROR: no valid entries found in ${CONFIG}" >&2
    exit 1
fi

echo "Config  : ${CONFIG}"
echo "Entries : ${NTOTAL} GWAS studies"
echo "──────────────────────────────────────────────────────────"

# ── Submit one job per dataset ────────────────────────────────────────────────
SUBMITTED=0
SKIPPED=0

for LINE in "${LINES[@]}"; do

    # Parse semicolon-delimited fields
    IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS MEM TIME \
        <<< "${LINE}"

    # Basic sanity check: all required fields must be non-empty
    if [[ -z "${INPUT_PATH}" || -z "${GWAS_NAME}" || -z "${POPULATION}" || \
          -z "${BUILD}"      || -z "${MEM}"        || -z "${TIME}" ]]; then
        echo "SKIP (malformed line — missing required field): ${LINE}" >&2
        (( SKIPPED++ )) || true
        continue
    fi

    JOB_ID=$(sbatch \
        --job-name="gwaslab_${GWAS_NAME}" \
        --nodes="${NODES}" --cpus-per-task="${CPUS}" \
        --mem="${MEM}" \
        --time="${TIME}" \
        --mail-type="${MAIL_TYPE}" --mail-user="${EMAIL}" \
        --output="${LOG_BASE}/${GWAS_NAME}_%j.out" \
        --error="${LOG_BASE}/${GWAS_NAME}_%j.err" \
        "$@" \
        "${WORKER_SCRIPT}" "${LINE}" \
        | awk '{print $NF}')

    echo "Submitted ${GWAS_NAME}  (mem=${MEM}, time=${TIME}) → job ${JOB_ID}"
    (( SUBMITTED++ )) || true

done

echo "──────────────────────────────────────────────────────────"
echo "Done: ${SUBMITTED} job(s) submitted, ${SKIPPED} line(s) skipped."

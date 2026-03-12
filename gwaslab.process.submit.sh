#!/usr/bin/env bash
#
# submit_gwaslab.sh — count valid entries in a config file and submit the
#                     gwaslab_array.sh SLURM array job.
#
# Usage:
#   bash submit_gwaslab.sh gwas_list.txt
#   bash submit_gwaslab.sh gwas_list.txt --partition=highmem --time=48:00:00
#
# Any extra arguments are forwarded directly to sbatch, allowing you to
# override directives from the command line (e.g. --mem, --time, --partition).
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARRAY_SCRIPT="${SCRIPT_DIR}/gwaslab_array.sh"

# ── Arguments ─────────────────────────────────────────────────────────────────
CONFIG="${1:?Usage: bash submit_gwaslab.sh <config.txt> [extra sbatch args]}"
shift   # remaining args forwarded to sbatch

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config file not found: ${CONFIG}" >&2
    exit 1
fi

if [[ ! -f "${ARRAY_SCRIPT}" ]]; then
    echo "ERROR: array script not found: ${ARRAY_SCRIPT}" >&2
    exit 1
fi

# Count valid lines (non-blank, non-comment)
N=$(grep -cve '^\s*#' -e '^\s*$' "${CONFIG}" || true)

if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: no valid entries found in ${CONFIG}" >&2
    exit 1
fi

# ── Create log directory ───────────────────────────────────────────────────────
mkdir -p logs

# ── Submit ────────────────────────────────────────────────────────────────────
echo "Config      : ${CONFIG}"
echo "Entries     : ${N} GWAS studies"
echo "Array range : 1–${N}"
echo ""

sbatch --array="1-${N}" "$@" "${ARRAY_SCRIPT}" "${CONFIG}"

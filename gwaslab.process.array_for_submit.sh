#!/usr/bin/env bash
#
# gwaslab.process.array_for_submit.sh — SLURM array job for the gwaslab processing pipeline
#
# Submit with:  bash gwaslab.process.submit.sh gwas_list.txt
# Or directly:  sbatch --array=1-N gwaslab.process.array_for_submit.sh gwas_list.txt
#
# Config file format: tab-separated, one GWAS per line (see gwas_list.txt.example)
#   COL1  full path to input GWAS file
#   COL2  GWAS name               → --gwas / output subdirectory name
#   COL3  population              → --population  (EUR / EAS / SAS / AFR / AMR / META)
#   COL4  build                   → --build       (18 / 19 / 38)
#   COL5  N total          (or .) → --n
#   COL6  N cases          (or .) → --n-cases
#   COL7  N controls       (or .) → --n-controls
#
# If all three of COL5–COL7 are provided, --force-n is added automatically.
#
# ─────────────────────────────────────────────────────────────────────────────

# ── SLURM directives ─────────────────────────────────────────────────────────
#SBATCH --job-name=gwaslab_process       # Job name
#SBATCH --output=/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo/gwaslab_process%A_%a.out     # Standard output log file
#SBATCH --error=/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo/gwaslab_process%A_%a.err   # Error log
#SBATCH --ntasks=1              # Number of tasks (always 1 for a single Python process)
#SBATCH --cpus-per-task=8       # CPU cores — keep in sync with --threads below
#SBATCH --mem=128G               # Total RAM per job (64G or 128G)
#SBATCH --time=04:00:00              # Time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=s.w.vanderlaan[at]gmail[dot]com      # Where to send mail

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION — adjust these paths for your HPC environment
# ─────────────────────────────────────────────────────────────────────────────
PYTHON_SCRIPT="/hpc/local/Rocky8/dhl_ec/software/gwas2cojo/gwaslab.process.py"
REF_DIR="/hpc/dhl_ec/data/references/gwaslab"
OUT_BASE="/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo"
CONDA_ENV="gwas2cojo"

# ─────────────────────────────────────────────────────────────────────────────
# Activate conda environment
# ─────────────────────────────────────────────────────────────────────────────
# Uncomment the appropriate block for your HPC:

# --- Option A: conda ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# --- Option B: module + conda ---
# module load miniconda3
# conda activate "${CONDA_ENV}"

# ─────────────────────────────────────────────────────────────────────────────
# Load config file
# ─────────────────────────────────────────────────────────────────────────────
CONFIG="${1:-gwas_list.txt}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config file not found: ${CONFIG}" >&2
    exit 1
fi

# Collect all valid lines (skip blank lines and lines starting with #)
mapfile -t LINES < <(grep -v '^\s*#' "${CONFIG}" | grep -v '^\s*$')

NTOTAL="${#LINES[@]}"
if [[ "${SLURM_ARRAY_TASK_ID}" -lt 1 || "${SLURM_ARRAY_TASK_ID}" -gt "${NTOTAL}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (1–${NTOTAL})" >&2
    exit 1
fi

# SLURM_ARRAY_TASK_ID is 1-indexed; bash arrays are 0-indexed
LINE="${LINES[$((SLURM_ARRAY_TASK_ID - 1))]}"

if [[ -z "${LINE}" ]]; then
    echo "ERROR: empty line for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Parse columns
# ─────────────────────────────────────────────────────────────────────────────
IFS=$'\t' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS \
    <<< "${LINE}"

# Derive --directory and --input from the full path
GWAS_DIR=$(dirname  "${INPUT_PATH}")
INPUT_FILE=$(basename "${INPUT_PATH}")

# ─────────────────────────────────────────────────────────────────────────────
# Build command
# ─────────────────────────────────────────────────────────────────────────────
CMD=(
    python "${PYTHON_SCRIPT}"
    --gwas        "${GWAS_NAME}"
    --input       "${INPUT_FILE}"
    --directory   "${GWAS_DIR}"
    --ref         "${REF_DIR}"
    --output      "${OUT_BASE}/${GWAS_NAME}"
    --population  "${POPULATION}"
    --build       "${BUILD}"
    --liftover
    --figures
    --threads     8
    --dbsnp
    --qc
    --cojo --cojo-pos --cojo-id rsid
    --leads
    --fill-eaf
)

# ── Optional N arguments ──────────────────────────────────────────────────────
N_SET=0; NCASES_SET=0; NCONTROLS_SET=0

if [[ -n "${N:-}"        && "${N}"        != "." ]]; then CMD+=(--n         "${N}");        N_SET=1;         fi
if [[ -n "${N_CASES:-}"  && "${N_CASES}"  != "." ]]; then CMD+=(--n-cases   "${N_CASES}");  NCASES_SET=1;    fi
if [[ -n "${N_CONTROLS:-}" && "${N_CONTROLS}" != "." ]]; then CMD+=(--n-controls "${N_CONTROLS}"); NCONTROLS_SET=1; fi

# --force-n only when all three N values are explicitly provided
if [[ "${N_SET}" -eq 1 && "${NCASES_SET}" -eq 1 && "${NCONTROLS_SET}" -eq 1 ]]; then
    CMD+=(--force-n)
fi

# ─────────────────────────────────────────────────────────────────────────────
# Log and run
# ─────────────────────────────────────────────────────────────────────────────
echo "========================================================"
echo "Job         : ${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
echo "GWAS        : ${GWAS_NAME}"
echo "Input       : ${INPUT_PATH}"
echo "Output      : ${OUT_BASE}/${GWAS_NAME}"
echo "Population  : ${POPULATION}  |  Build: ${BUILD}"
echo "N / Ncases / Ncontrols: ${N:-'.'} / ${N_CASES:-'.'} / ${N_CONTROLS:-'.'}"
echo "Command     : ${CMD[*]}"
echo "========================================================"

"${CMD[@]}"

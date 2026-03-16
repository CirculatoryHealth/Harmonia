#!/usr/bin/env bash
#
# gwaslab.process.array_for_submit.sh — SLURM worker script for one GWAS dataset
#
# Do NOT submit this script directly.
# Use gwaslab.process.submit.sh, which calls sbatch once per dataset with the
# correct --mem, --time, --output, and --error set per job.
#
# The script expects at least one argument: a semicolon-delimited config line
# from gwas_list.txt (COL1–COL9 required; COL10–COL11 MEM_LIGHT/TIME_LIGHT are
# read by gwaslab.process.submit_staged.sh and ignored here).
#
# ─────────────────────────────────────────────────────────────────────────────

# ── SLURM directives (fixed across all jobs) ─────────────────────────────────
#SBATCH --ntasks=1              # One Python process per job
#SBATCH --cpus-per-task=8       # CPU cores — keep in sync with --threads below
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.w.vanderlaan[at]gmail[dot]com
# --job-name, --mem, --time, --output, --error are set by gwaslab.process.submit.sh

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
# --- Option A: conda ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# --- Option B: module + conda ---
# module load miniconda3
# conda activate "${CONDA_ENV}"

# ─────────────────────────────────────────────────────────────────────────────
# Parse the semicolon-delimited config line passed as $1
# Format: INPUT_PATH;GWAS_NAME;POPULATION;BUILD;N;N_CASES;N_CONTROLS;MEM;TIME[;MEM_LIGHT;TIME_LIGHT]
# COL10–COL11 (MEM_LIGHT/TIME_LIGHT) are consumed by submit_staged.sh; ignored here.
# ─────────────────────────────────────────────────────────────────────────────
LINE="${1:?ERROR: no config line provided. Submit via gwaslab.process.submit.sh or gwaslab.process.submit_staged.sh}"
STAGE="${2:-all}"   # pipeline stage; defaults to 'all' (full end-to-end run)

IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS MEM TIME _ _ \
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
    --stage       "${STAGE}"
)

# ── Optional N arguments ──────────────────────────────────────────────────────
N_SET=0; NCASES_SET=0; NCONTROLS_SET=0

if [[ -n "${N:-}"          && "${N}"          != "." ]]; then CMD+=(--n         "${N}");        N_SET=1;        fi
if [[ -n "${N_CASES:-}"    && "${N_CASES}"    != "." ]]; then CMD+=(--n-cases   "${N_CASES}");  NCASES_SET=1;   fi
if [[ -n "${N_CONTROLS:-}" && "${N_CONTROLS}" != "." ]]; then CMD+=(--n-controls "${N_CONTROLS}"); NCONTROLS_SET=1; fi

# --force-n only when all three N values are explicitly provided
if [[ "${N_SET}" -eq 1 && "${NCASES_SET}" -eq 1 && "${NCONTROLS_SET}" -eq 1 ]]; then
    CMD+=(--force-n)
fi

# ─────────────────────────────────────────────────────────────────────────────
# Log and run
# ─────────────────────────────────────────────────────────────────────────────
echo "========================================================"
echo "Job         : ${SLURM_JOB_ID:-local}"
echo "GWAS        : ${GWAS_NAME}"
echo "Input       : ${INPUT_PATH}"
echo "Output      : ${OUT_BASE}/${GWAS_NAME}"
echo "Population  : ${POPULATION}  |  Build: ${BUILD}"
echo "Stage       : ${STAGE}"
echo "Memory      : ${MEM}  |  Time limit: ${TIME}"
echo "N / Ncases / Ncontrols: ${N:-'.'} / ${N_CASES:-'.'} / ${N_CONTROLS:-'.'}"
echo "Command     : ${CMD[*]}"
echo "========================================================"

"${CMD[@]}"

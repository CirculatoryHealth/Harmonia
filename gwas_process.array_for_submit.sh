#!/usr/bin/env bash
#
# gwas_process.array_for_submit.sh — SLURM worker script for one GWAS dataset
#
# Do NOT submit this script directly.
# Use gwas_process.submit.sh, which calls sbatch once per dataset with the
# correct --mem, --time, --output, and --error set per job.
#
# The script expects at least one argument: a semicolon-delimited config line
# from gwas_list.txt (COL1–COL9 required; COL10–COL11 MEM_LIGHT/TIME_LIGHT are
# read by gwas_process.submit_staged.sh and ignored here).
#
# ─────────────────────────────────────────────────────────────────────────────

# ── SLURM directives (fixed across all jobs) ─────────────────────────────────
#SBATCH --ntasks=1              # One Python process per job
#SBATCH --cpus-per-task=8       # CPU cores — keep in sync with --threads below
# --job-name, --mem, --time, --output, --error, --mail-type, --mail-user
# are all set by the calling submit script (submit.sh or submit_staged.sh)
# via the sbatch command line.

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration — loaded from gwas2cojo.conf.
# Copy gwas2cojo.conf.example → gwas2cojo.conf and fill in your paths once.
#
# When running as a SLURM job the submit script exports GWAS2COJO_CONF with
# the absolute path to gwas2cojo.conf.  SLURM copies this script to its own
# spool directory before execution, so BASH_SOURCE[0] points there rather
# than to the original script location — making a BASH_SOURCE-relative lookup
# unreliable.  GWAS2COJO_CONF is therefore the preferred source.
# ─────────────────────────────────────────────────────────────────────────────
if [[ -n "${GWAS2COJO_CONF:-}" && -f "${GWAS2COJO_CONF}" ]]; then
    CONF="${GWAS2COJO_CONF}"          # path exported by the submit script
else
    # Fallback: look next to this script (works for direct/local invocation)
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONF="${SCRIPT_DIR}/gwas2cojo.conf"
fi
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: gwas2cojo.conf not found (tried: ${CONF})." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    echo "       When submitting via SLURM make sure the submit script exports GWAS2COJO_CONF." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL

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
LINE="${1:?ERROR: no config line provided. Submit via gwas_process.submit.sh or gwas_process.submit_staged.sh}"
STAGE="${2:-all}"   # pipeline stage; defaults to 'all' (full end-to-end run)

# When running as a SLURM array job, pass the task ID as --chrom so the Python
# script processes a single chromosome shard.  Task IDs: 1-22=autosomes,
# 23=X, 24=Y, 25=nonPAR, 26=MT.  If the chromosome is absent the Python script
# exits 0 gracefully, satisfying afterok for the next array stage.
CHROM_ARG=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CHROM_ARG=(--chrom "${SLURM_ARRAY_TASK_ID}")
fi

IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS MEM TIME _ _ EXTRA_FLAGS \
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
    --ldsc
    --leads
    --fill-eaf
    --stage       "${STAGE}"
    "${CHROM_ARG[@]}"
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

# ── Optional per-study extra flags (COL12 in gwas_list.txt) ──────────────────
# Use '.' as a no-op placeholder.  Multiple flags are space-separated, e.g.:
#   --keep-multiallelic
#   --keep-multiallelic --no-figures
if [[ -n "${EXTRA_FLAGS:-}" && "${EXTRA_FLAGS}" != "." ]]; then
    read -ra _extra <<< "${EXTRA_FLAGS}"
    CMD+=("${_extra[@]}")
fi

# ─────────────────────────────────────────────────────────────────────────────
# Log and run
# ─────────────────────────────────────────────────────────────────────────────
echo "========================================================"
echo "Job         : ${SLURM_JOB_ID:-local}${SLURM_ARRAY_TASK_ID:+ [array task ${SLURM_ARRAY_TASK_ID}]}"
echo "GWAS        : ${GWAS_NAME}"
echo "Input       : ${INPUT_PATH}"
echo "Output      : ${OUT_BASE}/${GWAS_NAME}"
echo "Population  : ${POPULATION}  |  Build: ${BUILD}"
echo "Stage       : ${STAGE}"
echo "Chromosome  : ${SLURM_ARRAY_TASK_ID:-'(whole-genome)'}"
echo "Memory      : ${MEM}  |  Time limit: ${TIME}"
echo "N / Ncases / Ncontrols: ${N:-'.'} / ${N_CASES:-'.'} / ${N_CONTROLS:-'.'}"
echo "Extra flags : ${EXTRA_FLAGS:-.}"
echo "Command     : ${CMD[*]}"
echo "========================================================"

"${CMD[@]}"

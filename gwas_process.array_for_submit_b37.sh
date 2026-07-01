#!/usr/bin/env bash
#
# gwas_process.array_for_submit_b37.sh — SLURM worker for GRCh37/hg19 output
#
# Identical to gwas_process.array_for_submit.sh except:
#   • --liftover is NOT passed (no hg19→hg38 forward liftover)
#   • --output-build 19 is passed instead (triggers hg38→hg19 reverse liftover
#     for BUILD=38 inputs; BUILD=19/37 inputs pass through unchanged)
#
# Do NOT submit this script directly.
# Use gwas_process.submit_staged_b37.sh.
#
# ─────────────────────────────────────────────────────────────────────────────

# ── SLURM directives (fixed across all jobs) ─────────────────────────────────
#SBATCH --ntasks=1              # One Python process per job
#SBATCH --cpus-per-task=8       # CPU cores — keep in sync with --threads below
# --job-name, --mem, --time, --output, --error, --mail-type, --mail-user
# are all set by the calling submit script (submit_staged_b37.sh)
# via the sbatch command line.

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration
# ─────────────────────────────────────────────────────────────────────────────
if [[ -n "${GWAS2COJO_CONF:-}" && -f "${GWAS2COJO_CONF}" ]]; then
    CONF="${GWAS2COJO_CONF}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONF="${SCRIPT_DIR}/gwas2cojo.conf"
fi
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: gwas2cojo.conf not found (tried: ${CONF})." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL

# ─────────────────────────────────────────────────────────────────────────────
# Activate conda environment
# ─────────────────────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# ─────────────────────────────────────────────────────────────────────────────
# Parse the semicolon-delimited config line passed as $1
# ─────────────────────────────────────────────────────────────────────────────
LINE="${1:?ERROR: no config line provided. Submit via gwas_process.submit_staged_b37.sh}"
STAGE="${2:-all}"

CHROM_ARG=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CHROM_ARG=(--chrom "${SLURM_ARRAY_TASK_ID}")
fi

IFS=';' read -r INPUT_PATH GWAS_NAME POPULATION BUILD N N_CASES N_CONTROLS MEM TIME _ _ EXTRA_FLAGS \
    <<< "${LINE}"

GWAS_DIR=$(dirname  "${INPUT_PATH}")
INPUT_FILE=$(basename "${INPUT_PATH}")

# ─────────────────────────────────────────────────────────────────────────────
# Build command — b37 variant:
#   • --output-build 19  instead of --liftover
#     hg38 inputs → reverse-lifted to GRCh37/hg19
#     hg19/hg37 inputs → passed through unchanged
# ─────────────────────────────────────────────────────────────────────────────
CMD=(
    python "${PYTHON_SCRIPT}"
    --gwas        "${GWAS_NAME}"
    --input       "${INPUT_FILE}"
    --directory   "${GWAS_DIR}"
    --ref         "${REF_DIR}"
    --output      "${OUT_BASE}/b37/${GWAS_NAME}"
    --population  "${POPULATION}"
    --build       "${BUILD}"
    --output-build 19
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

if [[ "${N_SET}" -eq 1 && "${NCASES_SET}" -eq 1 && "${NCONTROLS_SET}" -eq 1 ]]; then
    CMD+=(--force-n)
fi

# ── Optional per-study extra flags (COL12 in gwas_list_b37.txt) ──────────────
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
echo "Output      : ${OUT_BASE}/b37/${GWAS_NAME}"
echo "Population  : ${POPULATION}  |  Build: ${BUILD}  →  output: hg19/GRCh37"
echo "Stage       : ${STAGE}"
echo "Chromosome  : ${SLURM_ARRAY_TASK_ID:-'(whole-genome)'}"
echo "Memory      : ${MEM}  |  Time limit: ${TIME}"
echo "N / Ncases / Ncontrols: ${N:-'.'} / ${N_CASES:-'.'} / ${N_CONTROLS:-'.'}"
echo "Extra flags : ${EXTRA_FLAGS:-.}"
echo "Command     : ${CMD[*]}"
echo "========================================================"

"${CMD[@]}"

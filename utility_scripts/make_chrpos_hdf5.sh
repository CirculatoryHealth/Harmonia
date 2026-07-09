#!/usr/bin/env bash
#SBATCH --job-name=make_chrpos_hdf5
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=make_chrpos_hdf5_%j.out
#SBATCH --error=make_chrpos_hdf5_%j.err
#SBATCH --mail-type=END,FAIL
#
# make_chrpos_hdf5.sh — SLURM wrapper for make_chrpos_hdf5.py
#
# Converts a dbSNP VCF into per-chromosome HDF5 files for rsID→CHR:POS lookup.
# Run once per build. Writes output into REF_DIR (from gwas2cojo.conf).
#
# Submit:
#   sbatch make_chrpos_hdf5.sh [--build hg19|hg38|all] [--vcf PATH] [--email addr]
#
# Examples:
#   sbatch make_chrpos_hdf5.sh
#   sbatch make_chrpos_hdf5.sh --build hg38
#   sbatch make_chrpos_hdf5.sh --build all --email s.w.vanderlaan-2@umcutrecht.nl
#   # When the dbSNP VCF is GCF-named and auto-detection fails:
#   sbatch make_chrpos_hdf5.sh --build hg19 \
#       --vcf /hpc/dhl_ec/data/references/gwaslab/GCF_000001405.25.gz
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
BUILD="hg19"
EMAIL=""
VCF=""

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build)  BUILD="$2";  shift 2 ;;
        --email)  EMAIL="$2";  shift 2 ;;
        --vcf)    VCF="$2";    shift 2 ;;
        --) shift; break ;;
        -*) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
        *) break ;;
    esac
done

# Apply mail-user dynamically if provided
if [[ -n "${EMAIL}" && -n "${SLURM_JOB_ID:-}" ]]; then
    scontrol update JobId="${SLURM_JOB_ID}" MailUser="${EMAIL}"
fi

# ── Load gwas2cojo.conf ───────────────────────────────────────────────────────
# Resolution order:
#   1. GWAS2COJO_CONF env var (explicit override)
#   2. SLURM_SUBMIT_DIR — the directory where sbatch was called (SLURM sets this
#      automatically; gwas2cojo.conf lives alongside harmonia.py there)
#   3. BASH_SOURCE fallback — works for direct local invocation but fails when
#      SLURM copies the script to its spool directory
if [[ -n "${GWAS2COJO_CONF:-}" && -f "${GWAS2COJO_CONF}" ]]; then
    CONF="${GWAS2COJO_CONF}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/gwas2cojo.conf" ]]; then
    CONF="${SLURM_SUBMIT_DIR}/gwas2cojo.conf"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    CONF="${SCRIPT_DIR}/../gwas2cojo.conf"
fi
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: gwas2cojo.conf not found (tried: ${CONF})" >&2
    echo "       Submit from the gwas2cojo directory, or export GWAS2COJO_CONF=/path/to/gwas2cojo.conf" >&2
    exit 1
fi
source "${CONF}"
# Sets: PYTHON_SCRIPT  REF_DIR  OUT_BASE  CONDA_ENV  EMAIL

# Derive the installation root from PYTHON_SCRIPT so the Python helper script
# is always found next to harmonia.py, regardless of where sbatch was called.
ROOTDIR="$(dirname "${PYTHON_SCRIPT}")"

# ── Activate conda ────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "========================================================"
echo "Job          : ${SLURM_JOB_ID:-local}"
echo "Build        : ${BUILD}"
echo "REF_DIR      : ${REF_DIR}"
echo "Threads      : ${SLURM_CPUS_PER_TASK:-8}"
if [[ -n "${VCF}" ]]; then
echo "VCF override : ${VCF}"
fi
echo "========================================================"

VCF_ARG=()
if [[ -n "${VCF}" ]]; then
    VCF_ARG=(--vcf "${VCF}")
fi

python "${ROOTDIR}/utility_scripts/make_chrpos_hdf5.py" \
    --ref-dir  "${REF_DIR}"               \
    --build    "${BUILD}"                  \
    --threads  "${SLURM_CPUS_PER_TASK:-8}" \
    --complevel 3                          \
    "${VCF_ARG[@]}"

echo "Done."

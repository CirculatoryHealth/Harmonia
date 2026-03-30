#!/usr/bin/env bash
#SBATCH --job-name=neale_addvariantinfo
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=neale_addvariantinfo_%A_%a.out
#SBATCH --error=neale_addvariantinfo_%A_%a.err
#SBATCH --mail-type=FAIL
#
# fix_neale_ukb_addvariantinfo.array.sh
# SLURM array job — adds rsid, chr, bp, ref, alt, info from variants.tsv.bgz
# to Neale Lab UKB GWAS summary files.
#
# Submit:
#   sbatch fix_neale_ukb_addvariantinfo.array.sh [--data-dir <dir>] [--email <addr>]
#   sbatch fix_neale_ukb_addvariantinfo.array.sh --data-dir /hpc/ukbiobank/neale/data \
#       --email s.w.vanderlaan-2@umcutrecht.nl
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
DATA_DIR="/hpc/ukbiobank/neale/data"
EMAIL=""

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        --email)    EMAIL="$2";    shift 2 ;;
        --) shift; break ;;
        -*) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
        *) break ;;
    esac
done

# Apply mail-user dynamically (SBATCH directives are static, scontrol updates at runtime)
if [[ -n "${EMAIL}" && -n "${SLURM_JOB_ID:-}" ]]; then
    scontrol update JobId="${SLURM_JOB_ID}" MailUser="${EMAIL}"
fi

VARIANTS="${DATA_DIR}/../information/variants.tsv.bgz"
OUT_DIR="${DATA_DIR}"

# ── File list (one per array task) ───────────────────────────────────────────
FILES=(
    ""   # index 0 unused — SLURM arrays are 1-based
    "${DATA_DIR}/H8_BPV.gwas.imputed_v3.both_sexes.tsv.bgz"
    "${DATA_DIR}/H8_LABYRINTHITIS.gwas.imputed_v3.both_sexes.tsv.bgz"
    "${DATA_DIR}/H8_VERTIGO.gwas.imputed_v3.both_sexes.tsv.bgz"
    "${DATA_DIR}/H81.gwas.imputed_v3.both_sexes.tsv.bgz"
)

GWAS_FILE="${FILES[${SLURM_ARRAY_TASK_ID}]}"
stem=$(basename "${GWAS_FILE%.tsv.bgz}")
OUT="${OUT_DIR}/${stem}.addvariantinfo.tsv.gz"

echo "========================================================"
echo "Job          : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Input        : ${GWAS_FILE}"
echo "Output       : ${OUT}"
echo "Variants     : ${VARIANTS}"
echo "========================================================"

if [[ ! -f "${GWAS_FILE}" ]]; then
    echo "ERROR: input file not found: ${GWAS_FILE}" >&2; exit 1
fi
if [[ ! -f "${VARIANTS}" ]]; then
    echo "ERROR: variants file not found: ${VARIANTS}" >&2; exit 1
fi

zcat "${GWAS_FILE}" | awk -v variants="${VARIANTS}" '
BEGIN {
    OFS = "\t"
    cmd = "zcat " variants
    while ((cmd | getline line) > 0) {
        n = split(line, a, "\t")
        if (a[1] == "variant") next
        # 1=variant 2=chr 3=pos 4=ref 5=alt 6=rsid 10=info
        vinfo[a[1]] = a[6] "\t" a[2] "\t" a[3] "\t" a[4] "\t" a[5] "\t" a[10]
    }
    close(cmd)
}
NR == 1 { print "rsid\tchr\tbp\tref\talt\tinfo\t" $0; next }
{
    if ($1 in vinfo) print vinfo[$1] "\t" $0
    else             print ".\t.\t.\t.\t.\t.\t" $0
}
' | gzip > "${OUT}"

echo "Done: ${OUT}"

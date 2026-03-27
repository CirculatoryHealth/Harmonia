#!/bin/bash
#SBATCH --job-name=fix_bp_chrpos
#SBATCH --output=fix_bp_chrpos_%j.out
#SBATCH --error=fix_bp_chrpos_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

INDIR="_BloodPressure/_UKB_BP"

FILES=(
    "UKB-ICBPmeta750k_SBPsummaryResults.txt.gz"
    "UKB-ICBPmeta750k_DBPsummaryResults.txt.gz"
    "UKB-ICBPmeta750k_PPsummaryResults.txt.gz"
)

for f in "${FILES[@]}"; do
    IN="${INDIR}/${f}"
    OUT="${INDIR}/${f%.txt.gz}_chrpos.txt.gz"

    echo "[$(date)] Processing: ${IN}"

    zcat "${IN}" \
        | awk 'BEGIN{OFS="\t"} NR==1{print "CHR","BP",$0} NR>1{split($1,a,":"); print a[1],a[2],$0}' \
        | gzip > "${OUT}"

    echo "[$(date)] Written:    ${OUT}"
done

echo "[$(date)] All done."

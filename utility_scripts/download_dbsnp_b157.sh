#!/bin/bash
#SBATCH --job-name=download_dbsnp_b157
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --output=/hpc/dhl_ec/data/references/gwaslab/download_dbsnp_b157_%j.out
#SBATCH --error=/hpc/dhl_ec/data/references/gwaslab/download_dbsnp_b157_%j.err

set -euo pipefail
cd /hpc/dhl_ec/data/references/gwaslab/

echo "=== Starting dbSNP b157 download: $(date) ==="

# ── hg38 ──────────────────────────────────────────────────────────────────
echo "--- Downloading b157 hg38 ---"
wget --tries=10 --continue --read-timeout=120 \
     -O GCF_000001405.40.gz \
     "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz"

bgzip -t GCF_000001405.40.gz && echo "hg38 VCF integrity OK" || { echo "hg38 VCF FAILED bgzip check"; exit 1; }
tabix -p vcf GCF_000001405.40.gz
echo "hg38 tbi size: $(du -sh GCF_000001405.40.gz.tbi)"
bcftools index -n GCF_000001405.40.gz 2>&1 || echo "(bcftools count requires CSI; tbi built, use tabix)"

# ── hg19 ──────────────────────────────────────────────────────────────────
echo "--- Downloading b157 hg19 ---"
wget --tries=10 --continue --read-timeout=120 \
     -O GCF_000001405.25.gz \
     "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.25.gz"

bgzip -t GCF_000001405.25.gz && echo "hg19 VCF integrity OK" || { echo "hg19 VCF FAILED bgzip check"; exit 1; }
tabix -p vcf GCF_000001405.25.gz
echo "hg19 tbi size: $(du -sh GCF_000001405.25.gz.tbi)"

echo "=== All downloads complete: $(date) ==="
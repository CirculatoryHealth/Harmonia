#!/bin/bash
#SBATCH --job-name=dbsnp_b151_download
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --output=dbsnp_b151_download_%j.log

# Download dbSNP b151 VCF files (hg38 and hg19) and rebuild tabix indices.
# Mirrors the style of download_dbsnp_b157.sh in the reference folder.
#
# Usage: sbatch utility_scripts/download_dbsnp_b151.sh
#
# Set REF_DIR to the directory where your reference VCFs live.
REF_DIR="/hpc/dhl_ec/data/references/gwaslab"   # ← adjust if needed

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "[$(date)] Starting dbSNP b151 downloads — REF_DIR=${REF_DIR}"

# ── hg38 ──────────────────────────────────────────────────────────────────────
echo "[$(date)] Downloading 00-All.vcf.gz (b151 hg38)..."
wget --tries=10 --continue --read-timeout=120 \
     "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz"

echo "[$(date)] Verifying and indexing 00-All.vcf.gz..."
bgzip -t 00-All.vcf.gz && echo "OK" || { echo "ERROR: file truncated or corrupt"; exit 1; }
tabix -p vcf 00-All.vcf.gz
echo "[$(date)] hg38 done — index: $(du -sh 00-All.vcf.gz.tbi | cut -f1)"

# ── hg19 ──────────────────────────────────────────────────────────────────────
echo "[$(date)] Downloading 00-All.vcf.gz (b151 hg19)..."
wget --tries=10 --continue --read-timeout=120 \
     "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz" \
     -O 00-All.hg19.vcf.gz

echo "[$(date)] Verifying and indexing 00-All.hg19.vcf.gz..."
bgzip -t 00-All.hg19.vcf.gz && echo "OK" || { echo "ERROR: file truncated or corrupt"; exit 1; }
tabix -p vcf 00-All.hg19.vcf.gz
echo "[$(date)] hg19 done — index: $(du -sh 00-All.hg19.vcf.gz.tbi | cut -f1)"

echo "[$(date)] All done."

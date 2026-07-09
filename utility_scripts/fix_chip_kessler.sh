#!/usr/bin/env bash
#SBATCH --job-name=fix_chip_kessler
#SBATCH --output=fix_chip_kessler_%j.out
#SBATCH --error=fix_chip_kessler_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Kessler2022 CHIP — extract relevant columns from the harmonised GWAS Catalog file.
# standard_error is NA throughout (Firth regression via REGENIE); it is intentionally
# omitted from output — harmonia.py derives SE from ci_upper/ci_lower instead.
# rsid is also omitted; harmonia.py re-assigns rsIDs via --dbsnp.
# gsub(/\r/, "") strips CRLF from the source file (rsid is the last column and would
# otherwise be stored as "rsid\r" in the awk hash, causing $h["rsid"]→$0 to dump
# the entire record into the output).
#
# Source columns used (12 output columns):
#   name                    → SNPID  (chr:pos:EA:NEA format, always populated; renamed so
#                                      harmonia.py recognises it via the "snpid" alias)
#   chromosome              → chromosome
#   base_pair_location      → base_pair_location  (GRCh38/hg38 — harmonised)
#   effect_allele           → effect_allele
#   other_allele            → other_allele
#   effect_allele_frequency → effect_allele_frequency
#   odds_ratio              → odds_ratio
#   ci_upper                → ci_upper
#   ci_lower                → ci_lower
#   p_value                 → p_value
#   num_cases               → num_cases
#   num_controls            → num_controls

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_CHIP/CHIP_EUR_Kessler2022/harmonised"
f="GCST90165267.h.tsv.gz"
in="${INDIR}/${f}"
out="${INDIR}/${f/.h.tsv.gz/.parsed.txt.gz}"

if [[ ! -f "${in}" ]]; then
    echo "ERROR: input file not found: ${in}" >&2
    exit 1
fi

echo "→ ${f}"
zcat "${in}" | awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1;i<=NF;i++) { gsub(/\r/, "", $i); h[$i]=i }
    print "SNPID","chromosome","base_pair_location",\
          "effect_allele","other_allele",\
          "effect_allele_frequency",\
          "odds_ratio","ci_upper","ci_lower",\
          "p_value","num_cases","num_controls"
    next
}
{
    gsub(/\r$/, "")
    print $h["name"],$h["chromosome"],$h["base_pair_location"],\
          $h["effect_allele"],$h["other_allele"],\
          $h["effect_allele_frequency"],\
          $h["odds_ratio"],$h["ci_upper"],$h["ci_lower"],\
          $h["p_value"],$h["num_cases"],$h["num_controls"]
}' | gzip > "${out}"
echo "   Lines: $(zcat "${out}" | wc -l)"

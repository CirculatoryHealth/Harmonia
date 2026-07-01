#!/usr/bin/env bash
#SBATCH --job-name=fix_chip_kessler
#SBATCH --output=fix_chip_kessler_%j.out
#SBATCH --error=fix_chip_kessler_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Kessler2022 CHIP — extract relevant columns from the harmonised GWAS Catalog file.
# standard_error is NA throughout (Firth regression via REGENIE); gwas_process.py
# will derive SE from ci_upper/ci_lower via the existing CI→SE path in correct_columns().
#
# Source columns used:
#   name                    → SNPID  (chr:pos:EA:NEA format, always populated)
#   rsid                    → rsid
#   chromosome              → chromosome
#   base_pair_location      → base_pair_location  (GRCh38/hg38 — harmonised)
#   effect_allele           → effect_allele
#   other_allele            → other_allele
#   odds_ratio              → odds_ratio
#   standard_error          → standard_error  (NA for all rows; CI→SE used instead)
#   ci_upper                → ci_upper
#   ci_lower                → ci_lower
#   effect_allele_frequency → effect_allele_frequency
#   p_value                 → p_value
#   num_cases               → num_cases
#   num_controls            → num_controls

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_CHIP/CHIP_EUR_Kessler2022/harmonised"

for f in GCST90165267.h.tsv.gz; do
    in="${INDIR}/${f}"
    out="${INDIR}/${f/.h.tsv.gz/.parsed.txt.gz}"

    if [[ ! -f "${in}" ]]; then
        echo "ERROR: input file not found: ${in}" >&2
        exit 1
    fi

    echo "→ ${f}"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"}
    NR==1 {
        for(i=1;i<=NF;i++) h[$i]=i
        print "SNPID","rsid","chromosome","base_pair_location",\
              "effect_allele","other_allele",\
              "odds_ratio","standard_error","ci_upper","ci_lower",\
              "effect_allele_frequency","p_value",\
              "num_cases","num_controls"
        next
    }
    {
        print $h["name"],$h["rsid"],$h["chromosome"],$h["base_pair_location"],\
              $h["effect_allele"],$h["other_allele"],\
              $h["odds_ratio"],$h["standard_error"],$h["ci_upper"],$h["ci_lower"],\
              $h["effect_allele_frequency"],$h["p_value"],\
              $h["num_cases"],$h["num_controls"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
done

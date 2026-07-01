#!/usr/bin/env bash
#SBATCH --job-name=fix_osa_verma
#SBATCH --output=fix_osa_verma_%j.out
#SBATCH --error=fix_osa_verma_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Verma2024 OSA — extract relevant columns; gwas_process.py derives SE from ci_upper/ci_lower
INDIR="/hpc/dhl_ec/data/_gwas_datasets/_OSA/Verma2024"

for f in GCST90475825.tsv.gz GCST90479998.tsv.gz; do
    in="${INDIR}/${f}"
    out="${INDIR}/${f/.tsv.gz/.parsed.txt.gz}"
    echo "→ ${f}"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"}
    NR==1 {
        for(i=1;i<=NF;i++) h[$i]=i
        print "chromosome","base_pair_location","effect_allele","other_allele",\
              "odds_ratio","ci_upper","ci_lower","effect_allele_frequency","p_value",\
              "rsid","n","num_cases","num_controls"
        next
    }
    {
        print $h["chromosome"],$h["base_pair_location"],\
              $h["effect_allele"],$h["other_allele"],\
              $h["odds_ratio"],$h["ci_upper"],$h["ci_lower"],$h["effect_allele_frequency"],\
              $h["p_value"],$h["rsid"],$h["n"],$h["num_cases"],$h["num_controls"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
done
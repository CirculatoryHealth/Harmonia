#!/usr/bin/env bash
#SBATCH --job-name=fix_prca_wang
#SBATCH --output=fix_prca_wang_%j.out
#SBATCH --error=fix_prca_wang_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Wang2023 PrCa — extract relevant columns


INDIR="/hpc/dhl_ec/data/_gwas_datasets/_PrCa/Wang2023"

for f in GCST90274713.h.tsv.gz GCST90274714.h.tsv.gz; do
    in="${INDIR}/${f}"
    out="${INDIR}/${f/.h.tsv.gz/.parsed.txt.gz}"
    echo "→ ${f}"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"}
    NR==1 {
        for(i=1;i<=NF;i++) h[$i]=i
        print "chromosome","base_pair_location","effect_allele","other_allele",\
              "beta","standard_error","effect_allele_frequency","p_value","rsid"
        next
    }
    {
        print $h["chromosome"],$h["base_pair_location"],\
              $h["effect_allele"],$h["other_allele"],\
              $h["beta"],$h["standard_error"],$h["effect_allele_frequency"],\
              $h["p_value"],$h["rsid"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
done
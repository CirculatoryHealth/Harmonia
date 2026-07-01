#!/usr/bin/env bash
#SBATCH --job-name=fix_t1d_mcgrail
#SBATCH --output=fix_t1d_mcgrail_%j.out
#SBATCH --error=fix_t1d_mcgrail_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# McGrail2026 T1D — extract relevant columns from the harmonised GWAS Catalog file.
# Original has 13 columns; keeping 10 key columns to avoid OOM during preprocess.
# hm_coordinate_conversion, hm_code are dropped (harmonisation metadata only).
#
# Source columns used:
#   chromosome            → chromosome   (alias → CHR ✓)
#   base_pair_location    → base_pair_location (alias → POS ✓; GRCh38/hg38 harmonised)
#   effect_allele         → effect_allele (alias → EA ✓)
#   other_allele          → other_allele  (alias → NEA ✓)
#   beta                  → beta          (alias → BETA ✓)
#   standard_error        → standard_error (alias → SE ✓)
#   effect_allele_frequency → effect_allele_frequency (alias → EAF ✓)
#   p_value               → p_value       (alias → P ✓)
#   rsid                  → rsid          (alias → rsID ✓)
#   n                     → n             (alias → N ✓)
#   variant_id            → variant_id    (alias → SNPID ✓; format: chr_pos_ref_alt)

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_DIABETES/T1D/mcgrail_2026/harmonised"

f="GCST90824163.h.tsv.gz"
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
    print "chromosome","base_pair_location","effect_allele","other_allele",\
          "beta","standard_error","effect_allele_frequency","p_value",\
          "rsid","n","variant_id"
    next
}
{
    print $h["chromosome"],$h["base_pair_location"],\
          $h["effect_allele"],$h["other_allele"],\
          $h["beta"],$h["standard_error"],$h["effect_allele_frequency"],\
          $h["p_value"],$h["rsid"],$h["n"],$h["variant_id"]
}' | gzip > "${out}"
echo "   Lines: $(zcat "${out}" | wc -l)"

#!/usr/bin/env bash
#SBATCH --job-name=fix_pgc_ptsd_nievergelt
#SBATCH --output=fix_pgc_ptsd_nievergelt_%j.out
#SBATCH --error=fix_pgc_ptsd_nievergelt_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Nievergelt2024 PTSD — convert PGC VCF format to a standard TSV.
# The source files are VCF-format with ## metadata headers and #CHROM header line.
# Columns: #CHROM ID POS A1 A2 FREQ NEFF Z P DIRE
# Z-score only (no BETA/SE); harmonia.py derives BETA/SE at merge stage via
# the Z + EAF + N formula (run_merge): SE = 1/sqrt(2*EAF*(1-EAF)*N), BETA = Z*SE.
#
# Column aliases already covered in harmonia.py (no changes needed):
#   CHROM → CHR    ID → SNPID    POS    A1 → EA    A2 → NEA
#   FREQ → EAF    NEFF → N    Z    P
# DIRE is dropped (per-cohort direction string, not needed downstream).

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_PGC/PTSD/Nievergelt2024"

for f in \
    eur_ptsd_pcs_v4_aug3_2021.vcf.gz \
    trans_ptsd_pcs_v4_aug3_2021.vcf.gz
do
    in="${INDIR}/${f}"
    out="${INDIR}/${f/.vcf.gz/.parsed.txt.gz}"

    if [[ ! -f "${in}" ]]; then
        echo "SKIP: file not found: ${in}" >&2
        continue
    fi

    echo "→ ${f}"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"; header=0}
    /^##/ { next }
    header==0 {
        $1 = ($1 ~ /^#/) ? substr($1, 2) : $1
        for(i=1;i<=NF;i++) h[$i]=i
        print "CHROM","ID","POS","A1","A2","FREQ","NEFF","Z","P"
        header=1; next
    }
    {
        print $h["CHROM"],$h["ID"],$h["POS"],$h["A1"],$h["A2"],\
              $h["FREQ"],$h["NEFF"],$h["Z"],$h["P"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
done

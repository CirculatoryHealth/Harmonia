#!/usr/bin/env bash
#SBATCH --job-name=fix_pgc_an_watson
#SBATCH --output=fix_pgc_an_watson_%j.out
#SBATCH --error=fix_pgc_an_watson_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Watson2019 AN — convert PGC VCF-TSV format to a standard TSV.
# File: pgcAN2.2019-07.vcf.tsv.gz (may or may not have ## VCF metadata headers).
# Columns: CHROM POS ID REF ALT BETA SE PVAL NGT IMPINFO NEFFDIV2 NCAS NCON DIRE
#
# Note: column order differs from PTSD/SCZ (POS before ID, REF/ALT instead of A1/A2).
# REF → NEA (alias "ref" ✓); ALT → EA (alias "alt" ✓)
# PVAL → PVAL (alias "pval" → P ✓); IMPINFO → IMPINFO (alias → INFO ✓)
# NCAS → NCAS (alias "ncas" → N_cases ✓); NCON → NCON (alias "ncon" → N_controls ✓)
# NEFFDIV2 → NEFF in output (multiply by 2 for effective N; alias "neff" → N ✓)
# NGT, DIRE dropped.
# N is NOT set in gwas_list.txt (dots) — NCAS/NCON columns provide per-variant N.

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_PGC/AN/Watson2019"

f="pgcAN2.2019-07.vcf.tsv.gz"
in="${INDIR}/${f}"
out="${INDIR}/${f/.vcf.tsv.gz/.parsed.txt.gz}"

if [[ ! -f "${in}" ]]; then
    echo "ERROR: input file not found: ${in}" >&2
    exit 1
fi

echo "→ ${f}"
zcat "${in}" | awk 'BEGIN{FS=OFS="\t"; header=0}
/^##/ { next }
header==0 {
    $1 = ($1 ~ /^#/) ? substr($1, 2) : $1
    for(i=1;i<=NF;i++) h[$i]=i
    print "CHROM","ID","POS","REF","ALT","BETA","SE","PVAL","IMPINFO","NEFF","NCAS","NCON"
    header=1; next
}
{
    neff = ($h["NEFFDIV2"] != "NA" && $h["NEFFDIV2"] != "") ? $h["NEFFDIV2"] * 2 : "NA"
    print $h["CHROM"],$h["ID"],$h["POS"],$h["REF"],$h["ALT"],\
          $h["BETA"],$h["SE"],$h["PVAL"],$h["IMPINFO"],\
          neff,$h["NCAS"],$h["NCON"]
}' | gzip > "${out}"
echo "   Lines: $(zcat "${out}" | wc -l)"

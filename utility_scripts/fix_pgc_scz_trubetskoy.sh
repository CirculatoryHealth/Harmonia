#!/usr/bin/env bash
#SBATCH --job-name=fix_pgc_scz_trubetskoy
#SBATCH --output=fix_pgc_scz_trubetskoy_%j.out
#SBATCH --error=fix_pgc_scz_trubetskoy_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Trubetskoy2022 SCZ — convert PGC VCF-TSV format to a standard TSV.
# The source files use VCF-style column names (may have ## metadata headers).
#
# EUR columns: CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL NCAS NCON NEFF
# PAN columns: CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL NGT DIRE NCAS NCON NEFFDIV2
#
# FCAS = A1 freq in cases; FCON = A1 freq in controls → use FCON as EAF proxy
# (renamed to EAF in output since "fcon" is not in harmonia.py EAF aliases).
# PVAL → PVAL (alias "pval" → P ✓); IMPINFO → IMPINFO (alias → INFO ✓)
# NCAS → NCAS (alias "ncas" → N_cases ✓); NCON → NCON (alias "ncon" → N_controls ✓)
# PAN: NEFFDIV2 → NEFF in output (multiply by 2 for effective N)
# EUR: NEFF kept as-is.
# NGT, DIRE, FCAS dropped.

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_PGC/SCZ/Trubetskoy2022"

# ── EUR ─────────────────────────────────────────────────────────────────────
f="PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
in="${INDIR}/${f}"
out="${INDIR}/${f/.vcf.tsv.gz/.parsed.txt.gz}"
if [[ -f "${in}" ]]; then
    echo "→ ${f} (EUR)"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"; header=0}
    /^##/ { next }
    header==0 {
        $1 = ($1 ~ /^#/) ? substr($1, 2) : $1
        for(i=1;i<=NF;i++) h[$i]=i
        print "CHROM","ID","POS","A1","A2","EAF","IMPINFO","BETA","SE","PVAL","NEFF","NCAS","NCON"
        header=1; next
    }
    {
        print $h["CHROM"],$h["ID"],$h["POS"],$h["A1"],$h["A2"],\
              $h["FCON"],$h["IMPINFO"],$h["BETA"],$h["SE"],$h["PVAL"],\
              $h["NEFF"],$h["NCAS"],$h["NCON"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
else
    echo "SKIP: file not found: ${in}" >&2
fi

# ── PAN ─────────────────────────────────────────────────────────────────────
f="PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz"
in="${INDIR}/${f}"
out="${INDIR}/${f/.vcf.tsv.gz/.parsed.txt.gz}"
if [[ -f "${in}" ]]; then
    echo "→ ${f} (PAN)"
    zcat "${in}" | awk 'BEGIN{FS=OFS="\t"; header=0}
    /^##/ { next }
    header==0 {
        $1 = ($1 ~ /^#/) ? substr($1, 2) : $1
        for(i=1;i<=NF;i++) h[$i]=i
        print "CHROM","ID","POS","A1","A2","EAF","IMPINFO","BETA","SE","PVAL","NEFF","NCAS","NCON"
        header=1; next
    }
    {
        # PAN uses NEFFDIV2; multiply by 2 to get effective N
        neff = ($h["NEFFDIV2"] != "NA" && $h["NEFFDIV2"] != "") ? $h["NEFFDIV2"] * 2 : "NA"
        print $h["CHROM"],$h["ID"],$h["POS"],$h["A1"],$h["A2"],\
              $h["FCON"],$h["IMPINFO"],$h["BETA"],$h["SE"],$h["PVAL"],\
              neff,$h["NCAS"],$h["NCON"]
    }' | gzip > "${out}"
    echo "   Lines: $(zcat "${out}" | wc -l)"
else
    echo "SKIP: file not found: ${in}" >&2
fi

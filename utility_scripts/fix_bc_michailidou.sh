#!/usr/bin/env bash
#SBATCH --job-name=fix_bc_michailidou
#SBATCH --output=fix_bc_michailidou_%j.out
#SBATCH --error=fix_bc_michailidou_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

# Michailidou2017 BC — split combined BCAC file into 3 separate analyses:
#   bc_all   : combined (OncoArray + iCOGS, all ancestry)
#   bc_erpos : ER-positive cases
#   bc_erneg : ER-negative cases
#
# Source columns used (same NEA/EA/CHR/POS for all three analyses):
#   var_name                              → SNPID
#   phase3_1kg_id                         → rsID
#   chr                                   → CHR
#   position_b37                          → POS
#   a0                                    → NEA (non-effect / reference allele)
#   a1                                    → EA  (effect / risk allele)
#   bcac_onco_icogs_gwas_eaf_controls     → EAF (combined)
#   bcac_onco_icogs_gwas_beta             → BETA (combined)
#   bcac_onco_icogs_gwas_se               → SE (combined)
#   bcac_onco_icogs_gwas_P1df             → P (combined)
#   bcac_onco_icogs_gwas_erpos_eaf_controls → EAF (ER+)
#   bcac_onco_icogs_gwas_erpos_beta         → BETA (ER+)
#   bcac_onco_icogs_gwas_erpos_se           → SE (ER+)
#   bcac_onco_icogs_gwas_erpos_P1df         → P (ER+)
#   bcac_onco_icogs_gwas_erneg_eaf_controls → EAF (ER-)
#   bcac_onco_icogs_gwas_erneg_beta         → BETA (ER-)
#   bcac_onco_icogs_gwas_erneg_se           → SE (ER-)
#   bcac_onco_icogs_gwas_erneg_P1df         → P (ER-)

set -euo pipefail

INDIR="/hpc/dhl_ec/data/_gwas_datasets/_BC/Michailidou2017"
INPUT="${INDIR}/oncoarray_bcac_public_release_oct17.txt.gz"

if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: input file not found: ${INPUT}" >&2
    exit 1
fi

echo "Input : ${INPUT}"
echo "──────────────────────────────────────────────────────────────"

# ── combined (all-BC) ────────────────────────────────────────────
OUT="${INDIR}/oncoarray_bcac_public_release_oct17.bc_all.txt.gz"
echo "→ bc_all"
zcat "${INPUT}" | awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1;i<=NF;i++) h[$i]=i
    print "SNPID","rsID","CHR","POS","NEA","EA","EAF","BETA","SE","P"
    next
}
{
    print $h["var_name"],$h["phase3_1kg_id"],$h["chr"],$h["position_b37"],\
          $h["a0"],$h["a1"],\
          $h["bcac_onco_icogs_gwas_eaf_controls"],\
          $h["bcac_onco_icogs_gwas_beta"],\
          $h["bcac_onco_icogs_gwas_se"],\
          $h["bcac_onco_icogs_gwas_P1df"]
}' | gzip > "${OUT}"
echo "   Lines: $(zcat "${OUT}" | wc -l)"

# ── ER-positive ───────────────────────────────────────────────────
OUT="${INDIR}/oncoarray_bcac_public_release_oct17.bc_erpos.txt.gz"
echo "→ bc_erpos"
zcat "${INPUT}" | awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1;i<=NF;i++) h[$i]=i
    print "SNPID","rsID","CHR","POS","NEA","EA","EAF","BETA","SE","P"
    next
}
{
    print $h["var_name"],$h["phase3_1kg_id"],$h["chr"],$h["position_b37"],\
          $h["a0"],$h["a1"],\
          $h["bcac_onco_icogs_gwas_erpos_eaf_controls"],\
          $h["bcac_onco_icogs_gwas_erpos_beta"],\
          $h["bcac_onco_icogs_gwas_erpos_se"],\
          $h["bcac_onco_icogs_gwas_erpos_P1df"]
}' | gzip > "${OUT}"
echo "   Lines: $(zcat "${OUT}" | wc -l)"

# ── ER-negative ───────────────────────────────────────────────────
OUT="${INDIR}/oncoarray_bcac_public_release_oct17.bc_erneg.txt.gz"
echo "→ bc_erneg"
zcat "${INPUT}" | awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1;i<=NF;i++) h[$i]=i
    print "SNPID","rsID","CHR","POS","NEA","EA","EAF","BETA","SE","P"
    next
}
{
    print $h["var_name"],$h["phase3_1kg_id"],$h["chr"],$h["position_b37"],\
          $h["a0"],$h["a1"],\
          $h["bcac_onco_icogs_gwas_erneg_eaf_controls"],\
          $h["bcac_onco_icogs_gwas_erneg_beta"],\
          $h["bcac_onco_icogs_gwas_erneg_se"],\
          $h["bcac_onco_icogs_gwas_erneg_P1df"]
}' | gzip > "${OUT}"
echo "   Lines: $(zcat "${OUT}" | wc -l)"

echo "──────────────────────────────────────────────────────────────"
echo "Done. Run fix_bc_michailidou_verify.sh (or zcat | head) to spot-check outputs."

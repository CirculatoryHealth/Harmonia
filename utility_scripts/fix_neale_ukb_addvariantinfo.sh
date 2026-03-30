#!/usr/bin/env bash
#
# fix_neale_ukb_addvariantinfo.sh
# Joins Neale Lab UKB GWAS summary files with the variants annotation file to
# add rsid, chr, bp, ref, alt, and info columns.
#
# Input  : {trait}.gwas.imputed_v3.both_sexes.tsv.bgz
# Output : <data-dir>/{trait}.gwas.imputed_v3.both_sexes.addvariantinfo.tsv.gz
#
# Usage:
#   bash fix_neale_ukb_addvariantinfo.sh --data-dir <dir> <file1.tsv.bgz> [file2.tsv.bgz ...]
#   bash fix_neale_ukb_addvariantinfo.sh --data-dir /ukbiobank/neale/data *.tsv.bgz
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
DATA_DIR="/hpc/dhl_ec/data/_gwas_datasets/_UKBB_Neale"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        --) shift; break ;;
        -*) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
        *) break ;;
    esac
done

VARIANTS="${DATA_DIR}/../information/variants.tsv.bgz"
OUT_DIR="${DATA_DIR}"

if [[ $# -eq 0 ]]; then
    echo "Usage: bash $(basename "$0") [--data-dir DIR] <file1.tsv.bgz> [file2.tsv.bgz ...]" >&2
    exit 1
fi

if [[ ! -f "${VARIANTS}" ]]; then
    echo "ERROR: variants file not found: ${VARIANTS}" >&2
    exit 1
fi

echo "Variants file : ${VARIANTS}"
echo "Output dir    : ${OUT_DIR}"
echo "Files to process: $#"
echo "──────────────────────────────────────────────────────"

for gwas_file in "$@"; do
    if [[ ! -f "${gwas_file}" ]]; then
        echo "SKIP (not found): ${gwas_file}" >&2
        continue
    fi

    stem=$(basename "${gwas_file%.tsv.bgz}")
    out="${OUT_DIR}/${stem}.addvariantinfo.tsv.gz"

    echo "Processing : ${stem}"

    zcat "${gwas_file}" | awk -v variants="${VARIANTS}" '
    BEGIN {
        OFS = "\t"
        # Load variant annotation: variant -> rsid, chr, pos, ref, alt, info
        # variants.tsv.bgz columns (1-based):
        #   1=variant  2=chr  3=pos  4=ref  5=alt  6=rsid  10=info
        cmd = "zcat " variants
        while ((cmd | getline line) > 0) {
            n = split(line, a, "\t")
            if (a[1] == "variant") continue       # skip header
            vinfo[a[1]] = a[6] "\t" a[2] "\t" a[3] "\t" a[4] "\t" a[5] "\t" a[10]
        }
        close(cmd)
    }
    NR == 1 {
        # Prepend new column headers before existing GWAS header
        print "rsid\tchr\tbp\tref\talt\tinfo\t" $0
        next
    }
    {
        if ($1 in vinfo)
            print vinfo[$1] "\t" $0
        else
            print ".\t.\t.\t.\t.\t.\t" $0
    }
    ' | gzip > "${out}"

    echo "  -> ${out}"
done

echo "──────────────────────────────────────────────────────"
echo "Done."

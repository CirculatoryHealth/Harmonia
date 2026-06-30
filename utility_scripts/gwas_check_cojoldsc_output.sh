#!/usr/bin/env bash
#
# gwas_check_cojoldsc_output.sh — check COJO and LDSC output files for gwas2cojo studies.
#
# For each study the script reports:
#   COJO_N : number of variants in the *.cojo (or *.cojo.gz) file,
#            MISS if the file is absent.
#   LDSC_N : number of variants in the *.ldsc.ldsc.tsv.gz file,
#            MISS if the file is absent.
#            Counts below 1,000   are flagged  ✗ (critically low).
#            Counts below 100,000 are flagged  ⚠  (suspicious).
#
# Usage:
#   gwas_check_outputs.sh [OPTIONS] [STUDY ...]
#   gwas_check_outputs.sh [OPTIONS] --list FILE
#
# Options:
#   -b, --base  DIR   Base directory that contains the per-study folders.
#                     Default: /hpc/dhl_ec/data/_gwas_datasets/gwas2cojo
#   -l, --list  FILE  Plain-text file of study names, one per line.
#                     Lines starting with # and blank lines are ignored.
#   -h, --help        Print this help and exit.
#
# Examples:
#   # Check a handful of studies in the default base directory
#   gwas_check_outputs.sh ADHD_EUR_Demontis2023 ALS_EUR_Rheenen2021
#
#   # Check all studies listed in a file, finished batch
#   gwas_check_outputs.sh --base /hpc/dhl_ec/data/_gwas_datasets/gwas2cojo/_finished_2026 \
#                         --list my_studies.txt
#

set -uo pipefail

# ── defaults ──────────────────────────────────────────────────────────────────
readonly SCRIPT_NAME="$(basename "$0")"
DEFAULT_BASE="/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo"
BASE="$DEFAULT_BASE"
LIST_FILE=""
STUDIES=()

# ── helpers ───────────────────────────────────────────────────────────────────
usage() {
    sed -n '/^# Usage/,/^[^#]/{ /^[^#]/d; s/^# \{0,1\}//; p }' "$0"
}

die() { echo "${SCRIPT_NAME}: error: $*" >&2; exit 1; }

# Count lines in a plain or gzip-compressed file, minus the header row.
count_variants() {
    local f="$1"
    local total
    if [[ "$f" == *.gz ]]; then
        total=$(zcat "$f" | wc -l)
    else
        total=$(wc -l < "$f")
    fi
    echo $(( total - 1 ))
}

# ── argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--base)
            [[ -n "${2:-}" ]] || die "--base requires an argument"
            BASE="$2"; shift 2 ;;
        -l|--list)
            [[ -n "${2:-}" ]] || die "--list requires an argument"
            LIST_FILE="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        -*)
            die "unknown option: $1" ;;
        *)
            STUDIES+=("$1"); shift ;;
    esac
done

# ── load study names from file if requested ────────────────────────────────────
if [[ -n "$LIST_FILE" ]]; then
    [[ -f "$LIST_FILE" ]] || die "list file not found: $LIST_FILE"
    while IFS= read -r line; do
        line="${line%%#*}"                    # strip inline comments
        line="${line#"${line%%[! ]*}"}"       # strip leading spaces
        line="${line%"${line##*[! ]}"}"       # strip trailing spaces
        [[ -n "$line" ]] && STUDIES+=("$line")
    done < "$LIST_FILE"
fi

[[ ${#STUDIES[@]} -gt 0 ]] || \
    die "no studies specified; use --list FILE or pass study names as arguments"

# ── output header ─────────────────────────────────────────────────────────────
FMT="%-55s  %12s  %12s  %s\n"
SEP_STUDY="$(printf '%0.s-' {1..55})"
SEP_N="$(printf '%0.s-' {1..12})"

printf "$FMT" "STUDY"      "COJO_N"  "LDSC_N"  "NOTE"
printf "$FMT" "$SEP_STUDY" "$SEP_N"  "$SEP_N"  "----"

# ── counters ──────────────────────────────────────────────────────────────────
n_total=0
n_cojo_miss=0
n_ldsc_ok=0
n_ldsc_low=0
n_ldsc_crit=0
n_ldsc_miss=0

# ── per-study check ───────────────────────────────────────────────────────────
for S in "${STUDIES[@]}"; do
    DIR="${BASE}/${S}"
    n_total=$(( n_total + 1 ))
    NOTE=""

    # ── COJO ──────────────────────────────────────────────────────────────────
    # Look for the main COJO data file: *.cojo or *.cojo.gz
    # (excludes GCTA auxiliary outputs like *.cojo.snplist, *.cojo.badsnps, etc.)
    COJO_FILE=$(ls "${DIR}"/*.cojo "${DIR}"/*.cojo.gz 2>/dev/null \
                | grep -E '\.cojo(\.gz)?$' | head -1)

    if [[ -n "$COJO_FILE" ]]; then
        COJO_N=$(count_variants "$COJO_FILE")
    else
        COJO_N="MISS"
        n_cojo_miss=$(( n_cojo_miss + 1 ))
    fi

    # ── LDSC ──────────────────────────────────────────────────────────────────
    LDSC_FILE=$(ls "${DIR}"/*.ldsc.tsv.gz 2>/dev/null | head -1)

    if [[ -n "$LDSC_FILE" ]]; then
        N_LDSC=$(count_variants "$LDSC_FILE")

        if   (( N_LDSC < 1000 )); then
            NOTE="✗ LDSC critically low"
            n_ldsc_crit=$(( n_ldsc_crit + 1 ))
        elif (( N_LDSC < 100000 )); then
            NOTE="⚠ LDSC low"
            n_ldsc_low=$(( n_ldsc_low + 1 ))
        else
            n_ldsc_ok=$(( n_ldsc_ok + 1 ))
        fi

        printf "$FMT" "$S" "$COJO_N" "$N_LDSC" "$NOTE"
    else
        n_ldsc_miss=$(( n_ldsc_miss + 1 ))
        printf "$FMT" "$S" "$COJO_N" "MISS" "$NOTE"
    fi

done

# ── summary footer ────────────────────────────────────────────────────────────
echo ""
printf "Studies : %d total\n" "$n_total"
printf "COJO    : %d missing\n" "$n_cojo_miss"
printf "LDSC    : %d OK (≥100k)  |  %d low (<100k)  |  %d critical (<1k)  |  %d missing\n" \
    "$n_ldsc_ok" "$n_ldsc_low" "$n_ldsc_crit" "$n_ldsc_miss"

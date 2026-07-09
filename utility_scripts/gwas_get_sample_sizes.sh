#!/usr/bin/env bash
#
# gwas_get_sample_sizes.sh — extract N, N_cases, N_controls from processed output.
#
# Reads from the first available file per study directory (in priority order):
#   *.gwaslab.qc.parquet  → fastest, has N / NCAS / NCON columns if present
#   *.gwaslab.parquet     → unfiltered version
#   *.gwaslab.qc.tsv.gz  → fallback
#   *.gwaslab.tsv.gz      → fallback
#
# For each study the script reports the maximum observed value of each column
# (most studies have a constant N per variant; max is robust to missing rows).
#
# Usage:
#   bash utility_scripts/gwas_get_sample_sizes.sh
#   bash utility_scripts/gwas_get_sample_sizes.sh --base /path/to/dir
#   bash utility_scripts/gwas_get_sample_sizes.sh --out my_sizes.tsv
#
# Output TSV columns:
#   study  N_total  N_cases  N_controls  source_file
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="${SCRIPT_DIR}/../harmonia.conf"
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: ${CONF} not found." >&2; exit 1
fi
# shellcheck source=../harmonia.conf.example
source "${CONF}"

BASE="${OUT_BASE}"
OUT="sample_sizes_$(date +%Y%m%d).tsv"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --base) BASE="${2:?--base requires a path}"; shift 2 ;;
        --out)  OUT="${2:?--out requires a path}";  shift 2 ;;
        --help|-h) grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

printf 'study\tN_total\tN_cases\tN_controls\tsource_file\n' > "${OUT}"

for study_dir in "${BASE}"/*/; do
    [[ -d "${study_dir}" ]] || continue
    study=$(basename "${study_dir}")

    # Find the best available output file.
    src=""
    for pattern in \
        "${study_dir}"*.gwaslab.qc.parquet \
        "${study_dir}"*.gwaslab.parquet \
        "${study_dir}"*.gwaslab.qc.tsv.gz \
        "${study_dir}"*.gwaslab.tsv.gz
    do
        for f in ${pattern}; do
            [[ -f "${f}" ]] && { src="${f}"; break 2; }
        done
    done

    if [[ -z "${src}" ]]; then
        printf '%s\t.\t.\t.\tnot_found\n' "${study}" >> "${OUT}"
        continue
    fi

    python3 - "${study}" "${src}" >> "${OUT}" <<'PYEOF'
import sys, os

study, src = sys.argv[1], sys.argv[2]

TARGET_COLS = ("N", "NCAS", "NCON")

def max_val(series):
    v = series.dropna()
    v = v[v > 0]
    return int(v.max()) if len(v) > 0 else "."

try:
    if src.endswith(".parquet"):
        import pyarrow.parquet as pq
        import pandas as pd
        avail = pq.read_schema(src).names
        cols  = [c for c in TARGET_COLS if c in avail]
        df    = pd.read_parquet(src, columns=cols) if cols else pd.DataFrame()
    else:
        import pandas as pd
        df = pd.read_csv(
            src, sep="\t",
            usecols=lambda c: c in TARGET_COLS,
            na_values=["NA", "."],
            low_memory=False,
        )

    n    = max_val(df["N"])    if "N"    in df.columns else "."
    ncas = max_val(df["NCAS"]) if "NCAS" in df.columns else "."
    ncon = max_val(df["NCON"]) if "NCON" in df.columns else "."
    print(f"{study}\t{n}\t{ncas}\t{ncon}\t{os.path.basename(src)}")

except Exception as exc:
    print(f"{study}\t.\t.\t.\tERROR: {exc}")
PYEOF
done

echo "Written: ${OUT}" >&2

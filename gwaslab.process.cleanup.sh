#!/usr/bin/env bash
#
# gwaslab.process.cleanup.sh — remove intermediate checkpoint files after a
#                               successful staged pipeline run.
#
# What is removed (per study output directory):
#   *.preprocess.parquet      handoff: preprocess → process-normalize
#   *.preprocess.json         metadata sidecar for the above
#   *.normalize.pkl           handoff: process-normalize → process-check-ref
#   *.checkref.pkl            handoff: process-check-ref → process-infer-strand
#   *.inferstrand.pkl         handoff: process-infer-strand → process-assign-rsid
#   *.assignrsid.pkl          handoff: process-assign-rsid → process-check-af
#   *.pkl      (raw, no .qc.) handoff: process-check-af → qc  (final raw pickle)
#
# What is NEVER removed:
#   *.qc.pkl                  needed by --stage cojo; also useful for quick reruns
#   *.parquet                 final output
#   *.tsv.gz                  final output
#   *.qc.parquet              final QC output
#   *.qc.tsv.gz               final QC output
#   *.cojo.gz / *.qc.cojo.gz  final COJO output
#   *.leads.tsv / *.qc.leads.tsv  lead-variant tables
#   *.log / *.gwaslab_process.log logs
#   PLOTS/                    all plot files
#
# Usage:
#   bash gwaslab.process.cleanup.sh --study CAD_Aragam
#   bash gwaslab.process.cleanup.sh --all
#   bash gwaslab.process.cleanup.sh --all --dry-run
#   bash gwaslab.process.cleanup.sh --study CAD_Aragam --keep-raw-pkl
#   bash gwaslab.process.cleanup.sh --config gwas_list.txt
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
OUT_BASE="/hpc/dhl_ec/data/_gwas_datasets/gwas2cojo"

# ─────────────────────────────────────────────────────────────────────────────
# Defaults
# ─────────────────────────────────────────────────────────────────────────────
DRY_RUN=0
KEEP_RAW_PKL=0     # set to 1 to preserve the final *.pkl (raw, non-QC) pickle
KEEP_QC_PKL=1      # always kept unless --remove-qc-pkl is passed
MODE=""            # "study", "all", or "config"
STUDY_NAME=""
CONFIG_FILE=""

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --study)         MODE="study"; STUDY_NAME="${2:?--study requires a name}"; shift 2 ;;
        --all)           MODE="all";   shift ;;
        --config)        MODE="config"; CONFIG_FILE="${2:?--config requires a path}"; shift 2 ;;
        --dry-run)       DRY_RUN=1;    shift ;;
        --keep-raw-pkl)  KEEP_RAW_PKL=1; shift ;;
        --remove-qc-pkl) KEEP_QC_PKL=0;  shift ;;
        --help|-h)       usage ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${MODE}" ]]; then
    echo "ERROR: specify --study <NAME>, --all, or --config <gwas_list.txt>" >&2
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Core cleanup function — operates on a single study output directory
# ─────────────────────────────────────────────────────────────────────────────
cleanup_study() {
    local study_dir="$1"
    local study_name="$2"

    if [[ ! -d "${study_dir}" ]]; then
        echo "  [${study_name}] Directory not found: ${study_dir} — skipping." >&2
        return
    fi

    # Build list of glob patterns to remove
    local patterns=(
        "${study_dir}/"*.preprocess.parquet
        "${study_dir}/"*.preprocess.json
        "${study_dir}/"*.normalize.pkl
        "${study_dir}/"*.checkref.pkl
        "${study_dir}/"*.inferstrand.pkl
        "${study_dir}/"*.assignrsid.pkl
    )

    # *.pkl that is NOT *.qc.pkl — the final raw pickle from process-check-af
    # We handle this separately because glob can't exclude a suffix directly.
    local raw_pkl_pattern="${study_dir}/"*.pkl

    local total_removed=0
    local total_size=0

    # ── Remove named intermediate checkpoints ─────────────────────────────────
    for pattern in "${patterns[@]}"; do
        for f in ${pattern}; do
            [[ -f "${f}" ]] || continue
            local sz
            sz=$(du -sh "${f}" 2>/dev/null | cut -f1)
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "  [DRY-RUN] would remove: ${f}  (${sz})"
            else
                echo "  [REMOVE] ${f}  (${sz})"
                rm -f "${f}"
            fi
            (( total_removed++ )) || true
        done
    done

    # ── Remove raw *.pkl (non-QC) unless --keep-raw-pkl ───────────────────────
    if [[ "${KEEP_RAW_PKL}" -eq 0 ]]; then
        for f in ${raw_pkl_pattern}; do
            [[ -f "${f}" ]] || continue
            # Skip *.qc.pkl files — those are the QC-filtered pickles
            [[ "${f}" == *.qc.pkl ]] && continue
            local sz
            sz=$(du -sh "${f}" 2>/dev/null | cut -f1)
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "  [DRY-RUN] would remove raw pickle: ${f}  (${sz})"
            else
                echo "  [REMOVE] raw pickle: ${f}  (${sz})"
                rm -f "${f}"
            fi
            (( total_removed++ )) || true
        done
    fi

    # ── Optionally remove *.qc.pkl ─────────────────────────────────────────────
    if [[ "${KEEP_QC_PKL}" -eq 0 ]]; then
        for f in "${study_dir}/"*.qc.pkl; do
            [[ -f "${f}" ]] || continue
            local sz
            sz=$(du -sh "${f}" 2>/dev/null | cut -f1)
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "  [DRY-RUN] would remove QC pickle: ${f}  (${sz})"
            else
                echo "  [REMOVE] QC pickle: ${f}  (${sz})"
                rm -f "${f}"
            fi
            (( total_removed++ )) || true
        done
    fi

    if [[ "${total_removed}" -eq 0 ]]; then
        echo "  [${study_name}] Nothing to remove (already clean)."
    else
        local verb="Removed"
        [[ "${DRY_RUN}" -eq 1 ]] && verb="Would remove"
        echo "  [${study_name}] ${verb} ${total_removed} file(s)."
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Mode: single study
# ─────────────────────────────────────────────────────────────────────────────
if [[ "${MODE}" == "study" ]]; then
    STUDY_DIR="${OUT_BASE}/${STUDY_NAME}/GWASCatalog"
    echo "Cleaning up: ${STUDY_NAME}"
    [[ "${DRY_RUN}" -eq 1 ]] && echo "(dry-run — nothing will be deleted)"
    echo "──────────────────────────────────────────────────────────"
    cleanup_study "${STUDY_DIR}" "${STUDY_NAME}"
    echo "──────────────────────────────────────────────────────────"
    echo "Done."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# Mode: all studies in OUT_BASE
# ─────────────────────────────────────────────────────────────────────────────
if [[ "${MODE}" == "all" ]]; then
    echo "Cleaning all studies under: ${OUT_BASE}"
    [[ "${DRY_RUN}" -eq 1 ]] && echo "(dry-run — nothing will be deleted)"
    echo "──────────────────────────────────────────────────────────"
    for study_dir in "${OUT_BASE}"/*/GWASCatalog; do
        [[ -d "${study_dir}" ]] || continue
        study_name=$(basename "$(dirname "${study_dir}")")
        echo "  Study: ${study_name}"
        cleanup_study "${study_dir}" "${study_name}"
    done
    echo "──────────────────────────────────────────────────────────"
    echo "Done."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# Mode: studies from config file
# ─────────────────────────────────────────────────────────────────────────────
if [[ "${MODE}" == "config" ]]; then
    if [[ ! -f "${CONFIG_FILE}" ]]; then
        echo "ERROR: config file not found: ${CONFIG_FILE}" >&2; exit 1
    fi

    mapfile -t LINES < <(grep -v '^\s*#' "${CONFIG_FILE}" | grep -v '^\s*$')
    NTOTAL="${#LINES[@]}"
    if [[ "${NTOTAL}" -eq 0 ]]; then
        echo "ERROR: no valid entries in ${CONFIG_FILE}" >&2; exit 1
    fi

    echo "Config     : ${CONFIG_FILE}  (${NTOTAL} studies)"
    [[ "${DRY_RUN}" -eq 1 ]] && echo "(dry-run — nothing will be deleted)"
    echo "──────────────────────────────────────────────────────────"

    for LINE in "${LINES[@]}"; do
        IFS=';' read -r INPUT_PATH GWAS_NAME _ <<< "${LINE}"
        [[ -z "${GWAS_NAME}" ]] && continue
        echo "  Study: ${GWAS_NAME}"
        cleanup_study "${OUT_BASE}/${GWAS_NAME}/GWASCatalog" "${GWAS_NAME}"
    done

    echo "──────────────────────────────────────────────────────────"
    echo "Done."
    exit 0
fi

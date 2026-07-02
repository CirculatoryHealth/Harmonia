#!/usr/bin/env bash
#
# gwas_process.cleanup.sh — remove intermediate checkpoint files after a
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
#   *.chr*.normalize.parquet  per-chromosome handoff: split → process-check-ref
#   *.chr*.checkref.parquet   per-chromosome handoff: process-check-ref → process-infer-strand
#   *.chr*.inferstrand.parquet per-chromosome handoff: process-infer-strand → process-assign-rsid
#   *.chr*.assignrsid.parquet per-chromosome handoff: process-assign-rsid → process-check-af
#   *.chr*.checkaf.parquet    per-chromosome handoff: process-check-af → merge
#
# What is NEVER removed:
#   *.parquet                 final output
#   *.tsv.gz                  final output
#   *.qc.parquet              final QC output
#   *.qc.tsv.gz               final QC output
#   *.cojo.gz / *.qc.cojo.gz  final COJO output
#   *.leads.tsv / *.qc.leads.tsv  lead-variant tables
#   *.log / *.gwas_process.log logs
#   PLOTS/                    all plot files
#
# Log archiving (default: enabled):
#   SLURM *.out / *.err files found in LOG_DIR (default: OUT_BASE) that belong to
#   the study being cleaned are moved into ${OUT_BASE}/<STUDY>/logs/.
#   The logs/ directory is then compressed to logs.tar.gz and removed.
#   This keeps the submit directory tidy and preserves logs for gwas_process.check.py.
#   Disable with --no-archive-logs.  Override the source dir with --log-dir PATH.
#
# Usage:
#   bash gwas_process.cleanup.sh --study CAD_Aragam
#   bash gwas_process.cleanup.sh --all
#   bash gwas_process.cleanup.sh --all --dry-run
#   bash gwas_process.cleanup.sh --study CAD_Aragam --keep-normalize-pkl
#   bash gwas_process.cleanup.sh --study CAD_Aragam --keep-raw-pkl
#   bash gwas_process.cleanup.sh --study CAD_Aragam --keep-qc-pkl
#   bash gwas_process.cleanup.sh --config gwas_list.txt
#   bash gwas_process.cleanup.sh --study CAD_Aragam --no-archive-logs
#   bash gwas_process.cleanup.sh --all --log-dir /path/to/slurm/logs
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Site configuration — loaded from gwas2cojo.conf (next to this script).
# Copy gwas2cojo.conf.example → gwas2cojo.conf and fill in your paths once.
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="${SCRIPT_DIR}/gwas2cojo.conf"
if [[ ! -f "${CONF}" ]]; then
    echo "ERROR: ${CONF} not found." >&2
    echo "       Copy gwas2cojo.conf.example to gwas2cojo.conf and fill in your paths." >&2
    exit 1
fi
# shellcheck source=gwas2cojo.conf.example
source "${CONF}"
# OUT_BASE is set from gwas2cojo.conf

# ─────────────────────────────────────────────────────────────────────────────
# Defaults
# ─────────────────────────────────────────────────────────────────────────────
DRY_RUN=0
KEEP_NORMALIZE_PKL=0  # set to 1 via --keep-normalize-pkl to preserve *.normalize.pkl
KEEP_RAW_PKL=0        # set to 1 via --keep-raw-pkl to preserve the final raw (non-QC) pickle
KEEP_QC_PKL=0         # set to 1 via --keep-qc-pkl to preserve *.qc.pkl
ARCHIVE_LOGS=1     # move SLURM *.out/*.err into <study_dir>/logs/ (disable: --no-archive-logs)
LOG_DIR=""         # source dir for log archiving; defaults to OUT_BASE after conf is sourced
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
        --dry-run)             DRY_RUN=1;            shift ;;
        --keep-normalize-pkl)  KEEP_NORMALIZE_PKL=1; shift ;;
        --keep-raw-pkl)        KEEP_RAW_PKL=1;       shift ;;
        --keep-qc-pkl)         KEEP_QC_PKL=1;        shift ;;
        --no-archive-logs)     ARCHIVE_LOGS=0;       shift ;;
        --log-dir)          LOG_DIR="${2:?--log-dir requires a path}"; shift 2 ;;
        --help|-h)          usage ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${MODE}" ]]; then
    echo "ERROR: specify --study <NAME>, --all, or --config <gwas_list.txt>" >&2
    exit 1
fi

# Default log source directory: same as OUT_BASE (where submit scripts write SLURM logs)
[[ -z "${LOG_DIR}" ]] && LOG_DIR="${OUT_BASE}"

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
        # Whole-genome checkpoints
        "${study_dir}/"*.preprocess.parquet
        "${study_dir}/"*.preprocess.json
        "${study_dir}/"*.checkref.pkl
        "${study_dir}/"*.inferstrand.pkl
        "${study_dir}/"*.assignrsid.pkl
        # Per-chromosome intermediate parquets (one per stage per chromosome)
        "${study_dir}/"*.chr*.normalize.parquet
        "${study_dir}/"*.chr*.checkref.parquet
        "${study_dir}/"*.chr*.inferstrand.parquet
        "${study_dir}/"*.chr*.assignrsid.parquet
        "${study_dir}/"*.chr*.checkaf.parquet
    )

    # *.pkl that is NOT *.qc.pkl — the final raw pickle from process-check-af
    # We handle this separately because glob can't exclude a suffix directly.
    local raw_pkl_pattern="${study_dir}/*.pkl"

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

    # ── Remove *.normalize.pkl unless --keep-normalize-pkl ────────────────────
    if [[ "${KEEP_NORMALIZE_PKL}" -eq 0 ]]; then
        for f in "${study_dir}/"*.normalize.pkl; do
            [[ -f "${f}" ]] || continue
            local sz
            sz=$(du -sh "${f}" 2>/dev/null | cut -f1)
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "  [DRY-RUN] would remove normalize pickle: ${f}  (${sz})"
            else
                echo "  [REMOVE] normalize pickle: ${f}  (${sz})"
                rm -f "${f}"
            fi
            (( total_removed++ )) || true
        done
    fi

    # ── Remove raw *.pkl (non-QC, non-normalize) unless --keep-raw-pkl ────────
    if [[ "${KEEP_RAW_PKL}" -eq 0 ]]; then
        for f in ${raw_pkl_pattern}; do
            [[ -f "${f}" ]] || continue
            # Skip *.qc.pkl and *.normalize.pkl — handled separately above
            [[ "${f}" == *.qc.pkl       ]] && continue
            [[ "${f}" == *.normalize.pkl ]] && continue
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

    # ── Remove *.qc.pkl unless --keep-qc-pkl ──────────────────────────────────
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

    # ── Archive SLURM *.out / *.err into <study_dir>/logs/ ────────────────────
    local logs_dir="${study_dir%/}/logs"
    if [[ "${ARCHIVE_LOGS}" -eq 1 ]]; then
        local n_archived=0
        for f in "${LOG_DIR}/${study_name}_"*.out "${LOG_DIR}/${study_name}_"*.err; do
            [[ -f "${f}" ]] || continue
            if [[ "${DRY_RUN}" -eq 1 ]]; then
                echo "  [DRY-RUN] would archive: $(basename "${f}") → ${logs_dir}/"
            else
                mkdir -p "${logs_dir}"
                mv "${f}" "${logs_dir}/"
                echo "  [ARCHIVE] $(basename "${f}") → ${logs_dir}/"
            fi
            (( n_archived++ )) || true
        done
        if [[ "${n_archived}" -gt 0 ]]; then
            local verb="Archived"
            [[ "${DRY_RUN}" -eq 1 ]] && verb="Would archive"
            echo "  [${study_name}] ${verb} ${n_archived} log file(s) → ${logs_dir}/"
            echo "  [${study_name}] Run check.py against archived logs:"
            echo "    python gwas_process.check.py ${study_name} ${logs_dir}"
        elif [[ "${DRY_RUN}" -eq 0 ]]; then
            echo "  [${study_name}] No SLURM log files found in ${LOG_DIR} — nothing archived."
        fi
    fi

    # ── Compress logs/ → logs.tar.gz ─────────────────────────────────────────
    if [[ -d "${logs_dir}" ]]; then
        if [[ "${DRY_RUN}" -eq 1 ]]; then
            echo "  [DRY-RUN] would compress: ${logs_dir}/ → ${logs_dir}.tar.gz"
        else
            tar -czf "${logs_dir}.tar.gz" -C "${study_dir%/}" logs \
                && rm -rf "${logs_dir}" \
                && echo "  [COMPRESS] logs/ → ${logs_dir}.tar.gz"
        fi
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Mode: single study
# ─────────────────────────────────────────────────────────────────────────────
if [[ "${MODE}" == "study" ]]; then
    STUDY_DIR="${OUT_BASE}/${STUDY_NAME}"
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
    for study_dir in "${OUT_BASE}"/*/; do
        [[ -d "${study_dir}" ]] || continue
        study_name=$(basename "${study_dir}")
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
    if [[ ! -e "${CONFIG_FILE}" ]]; then
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
        cleanup_study "${OUT_BASE}/${GWAS_NAME}" "${GWAS_NAME}"
    done

    echo "──────────────────────────────────────────────────────────"
    echo "Done."
    exit 0
fi

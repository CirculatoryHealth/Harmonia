#!/usr/bin/env bash
# fix_withmultiallelic_logs.sh
# Moves log files that contain _withmultiallelic_ in their name
# but are sitting in the base study directory, back to the correct one.
#
# Usage: bash fix_withmultiallelic_logs.sh <LOG_BASE> [--dry-run]

LOG_BASE="${1:?Usage: bash fix_withmultiallelic_logs.sh <LOG_BASE> [--dry-run]}"
DRY_RUN=0
[[ "${2:-}" == "--dry-run" ]] && DRY_RUN=1

moved=0; skipped=0; errors=0

for base_dir in "${LOG_BASE}"/*/; do
    study=$(basename "${base_dir%/}")
    # Only process non-withmultiallelic study directories
    [[ "${study}" == *_withmultiallelic* ]] && continue

    # Search inside <study>/logs/ if it exists, otherwise directly in <study>/
    if [[ -d "${base_dir}logs" ]]; then
        search_dir="${base_dir}logs"
    else
        search_dir="${base_dir%/}"
    fi

    # Find files whose names contain _withmultiallelic_
    while IFS= read -r -d '' file; do
        fname=$(basename "${file}")
        # Derive correct study name: strip everything after _withmultiallelic
        correct_study="${fname%%_withmultiallelic_*}_withmultiallelic"
        correct_dir="${LOG_BASE}/${correct_study}/logs"

        if [[ ! -d "${LOG_BASE}/${correct_study}" ]]; then
            echo "WARN: no study dir '${LOG_BASE}/${correct_study}' — skipping: ${fname}" >&2
            (( skipped++ )) || true
            continue
        fi

        if [[ "${DRY_RUN}" -eq 1 ]]; then
            echo "DRY-RUN: ${file}  ->  ${correct_dir}/"
        else
            mkdir -p "${correct_dir}"
            echo "MOVE: ${file}  ->  ${correct_dir}/"
            mv "${file}" "${correct_dir}/" || { echo "ERROR moving ${file}" >&2; (( errors++ )) || true; continue; }
        fi
        (( moved++ )) || true
    done < <(find "${search_dir}" -maxdepth 1 -name "*_withmultiallelic_*" -print0)
done

echo "────────────────────────────────────────────"
[[ "${DRY_RUN}" -eq 1 ]] && echo "DRY-RUN — no files moved." || true
echo "Moved: ${moved}  Skipped: ${skipped}  Errors: ${errors}"

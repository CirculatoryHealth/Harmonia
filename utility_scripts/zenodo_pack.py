#!/usr/bin/env python3
"""
zenodo_pack.py – Create per-study .tar.gz archives from Harmonia output for Zenodo.

For each redistributable study in SupTab1, this script finds the corresponding
Harmonia output directory under OUT_BASE, collects the final output files
(COJO, LDSC, harmonised TSV.GZ, and metadata JSON — not intermediate
checkpoints), and writes {ABBREVIATION}.tar.gz into the target directory.

Prerequisites
-------------
    pip install openpyxl

Usage
-----
    python utility_scripts/zenodo_pack.py \\
        --excel   /path/to/Tables.xlsx \\
        --out-base /path/to/harmonia/output \\
        --dest    /path/to/per_study_archives \\
        [--dry-run] \\
        [--abbr LDL_EUR HDL_EUR ...]        # limit to specific studies

What gets packed
----------------
    Included (final outputs):
        *.cojo.gz            COJO-format summary statistics
        *.ldsc.tsv.gz        LDSC-ready munged summary statistics
        *.tsv.gz             Harmonised full summary statistics
        *.json               Metadata / build info

    Excluded (intermediate checkpoints):
        *.pkl                gwaslab pickle checkpoints
        *.parquet            per-chromosome parquet splits
        *.png  *.pdf         QC plots (not required for downstream analysis)
        slurm-*.out          SLURM log files

Column indices in SupTab1 (header at row "Trait_Number"):
    2  Abbreviation
    17 Public_redistribution
    21 Zenodo_upload          (non-empty = redistributable)
"""

import argparse
import os
import sys
import tarfile
from pathlib import Path

import openpyxl

# ── file-selection rules ─────────────────────────────────────────────────────
INCLUDE_SUFFIXES = {
    '.cojo.gz', '.ldsc.tsv.gz', '.tsv.gz', '.json',
}
EXCLUDE_SUFFIXES = {
    '.pkl', '.parquet', '.png', '.pdf',
}
EXCLUDE_PREFIXES = ('slurm-',)


def should_include(path: Path) -> bool:
    name = path.name
    if any(name.startswith(p) for p in EXCLUDE_PREFIXES):
        return False
    # Check exclude suffixes first (some overlap with include, e.g. .tsv.gz vs .ldsc.tsv.gz)
    for suf in EXCLUDE_SUFFIXES:
        if name.endswith(suf):
            return False
    for suf in INCLUDE_SUFFIXES:
        if name.endswith(suf):
            return True
    return False


# ── Excel reader ─────────────────────────────────────────────────────────────

def load_redistributable_abbrs(excel_path: str) -> list[str]:
    """Return abbreviations of all studies marked as redistributable in SupTab1."""
    wb = openpyxl.load_workbook(excel_path)
    ws = wb['SupTab1']

    header_row = None
    for i, row in enumerate(ws.iter_rows(values_only=True), 1):
        if row[0] == 'Trait_Number':
            header_row = i
            break
    if header_row is None:
        sys.exit('ERROR: could not find header row (Trait_Number) in SupTab1')

    abbrs = []
    for row in ws.iter_rows(min_row=header_row + 1, values_only=True):
        if not row[0]:
            continue
        abbr   = row[2]
        pub    = str(row[17] or '')
        upload = row[21]
        if abbr and 'Free' in pub and upload:
            abbrs.append(abbr)
    return abbrs


# ── packing ──────────────────────────────────────────────────────────────────

def pack_study(abbr: str, out_base: Path, dest: Path, dry_run: bool) -> bool:
    """
    Pack one study's Harmonia output into {abbr}.tar.gz in dest/.
    Returns True on success, False if the study directory was not found.
    """
    study_dir = out_base / abbr
    if not study_dir.is_dir():
        print(f'  WARNING: {study_dir} not found — skipping {abbr}')
        return False

    files = sorted(f for f in study_dir.rglob('*') if f.is_file() and should_include(f))
    if not files:
        print(f'  WARNING: no packable files found in {study_dir} — skipping {abbr}')
        return False

    archive = dest / f'{abbr}.tar.gz'
    total_mb = sum(f.stat().st_size for f in files) / 1e6

    print(f'  {abbr}  ({len(files)} files, {total_mb:.1f} MB) → {archive.name}',
          flush=True)
    for f in files:
        print(f'      {f.relative_to(study_dir)}')

    if dry_run:
        return True

    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'w:gz') as tar:
        for f in files:
            tar.add(f, arcname=f.relative_to(study_dir))

    packed_mb = archive.stat().st_size / 1e6
    print(f'      → {packed_mb:.1f} MB compressed')
    return True


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--excel',    required=True,
                        help='Path to Tables.xlsx (SupTab1 sheet)')
    parser.add_argument('--out-base', required=True,
                        help='Harmonia base output directory (OUT_BASE in harmonia.conf)')
    parser.add_argument('--dest',     required=True,
                        help='Destination directory for the .tar.gz archives')
    parser.add_argument('--abbr',     nargs='+', metavar='ABBREVIATION',
                        help='Pack only these study abbreviations (default: all redistributable)')
    parser.add_argument('--dry-run',  action='store_true',
                        help='Show what would be packed without writing anything')
    args = parser.parse_args()

    out_base = Path(args.out_base)
    if not out_base.is_dir():
        sys.exit(f'ERROR: --out-base directory not found: {out_base}')

    dest = Path(args.dest)

    print(f'Reading SupTab1 from {args.excel} …')
    all_abbrs = load_redistributable_abbrs(args.excel)
    print(f'Found {len(all_abbrs)} redistributable studies.\n')

    if args.abbr:
        abbrs = [a for a in args.abbr if a in set(all_abbrs)]
        missing = set(args.abbr) - set(all_abbrs)
        if missing:
            print(f'WARNING: these abbreviations not found / not redistributable: {missing}')
    else:
        abbrs = all_abbrs

    ok = 0
    skipped = 0
    for abbr in abbrs:
        if pack_study(abbr, out_base, dest, dry_run=args.dry_run):
            ok += 1
        else:
            skipped += 1
        print()

    print(f'Done.  Packed: {ok}  Skipped: {skipped}')
    if args.dry_run:
        print('(dry-run — no files were written)')


if __name__ == '__main__':
    main()

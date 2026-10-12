#!/usr/bin/env python3
"""
zenodo_upload.py – Upload Harmonia GWAS summary-statistics archives to Zenodo.

Reads SupTab1 from Tables.xlsx to discover which files belong to which Zenodo
deposit (identified by the DOI in column U), then streams each .tar.gz archive
and the per-upload manifest TSV to the appropriate record via the Zenodo Files
API.  Already-uploaded files are skipped automatically (idempotent).

Prerequisites
-------------
    pip install openpyxl requests

Usage
-----
    export ZENODO_TOKEN=<your-personal-access-token>

    # Dry-run (no uploads, just show what would happen)
    python zenodo_upload.py --excel Tables.xlsx --per-study /path/to/per_study --dry-run

    # Upload everything not yet uploaded
    python zenodo_upload.py --excel Tables.xlsx --per-study /path/to/per_study

    # Upload specific folders only
    python zenodo_upload.py --excel Tables.xlsx --per-study /path/to/per_study \\
        --upload upload17 upload18

    # Upload and publish each record when all its files are done
    python zenodo_upload.py --excel Tables.xlsx --per-study /path/to/per_study --publish

Column indices in SupTab1 (header at row 17, data from row 18):
    0  Trait_Number          9  Consortium
    1  Trait_or_disease     10  PMID
    2  Abbreviation         11  Reference
    3  Category             17  Public_redistribution
    4  Ancestry             18  Zenodo_md5
    5  Build                19  Zenodo_title
    6  Sample_size          20  DOI              ← column U
    7  Cases                21  Zenodo_upload    ← column V
    8  Controls             22  Uploaded         ← column W
"""

import argparse
import os
import re
import sys
from pathlib import Path

import openpyxl
import requests
from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper

# ── constants ────────────────────────────────────────────────────────────────
ZENODO_API   = 'https://zenodo.org/api'
CHUNK_SIZE   = 8 * 1024 * 1024   # 8 MB streaming chunks

HEADER_ROW   = 17   # 1-indexed row that contains column names
COL_ABBR     = 2
COL_PUB      = 17
COL_DOI      = 20
COL_UPLOAD   = 21
COL_UPLOADED = 22


# ── helpers ──────────────────────────────────────────────────────────────────

def deposit_id_from_doi(doi: str) -> str:
    """Extract numeric Zenodo deposit ID from a DOI like 10.5281/zenodo.12345678."""
    m = re.search(r'zenodo\.(\d+)', doi)
    if not m:
        raise ValueError(f'Cannot parse deposit ID from DOI: {doi!r}')
    return m.group(1)


def get_bucket_url(dep_id: str, token: str) -> str:
    """Return the S3-style bucket URL for a deposit."""
    url  = f'{ZENODO_API}/deposit/depositions/{dep_id}'
    resp = requests.get(url, headers={'Authorization': f'Bearer {token}'})
    resp.raise_for_status()
    return resp.json()['links']['bucket']


def list_bucket_files(bucket_url: str, token: str) -> dict[str, int]:
    """Return {filename: size_in_bytes} for files already present in the deposit bucket."""
    resp = requests.get(bucket_url, headers={'Authorization': f'Bearer {token}'})
    resp.raise_for_status()
    return {f['key']: f['size'] for f in resp.json().get('contents', [])}


def upload_file(bucket_url: str, local_path: Path, token: str, dry_run: bool) -> None:
    """Stream a local file to the Zenodo bucket via PUT with a progress bar."""
    fname = local_path.name
    size  = local_path.stat().st_size
    print(f'    ↑ {fname}  ({size / 1e9:.2f} GB)', flush=True)
    if dry_run:
        return

    with tqdm(total=size, unit='B', unit_scale=True, unit_divisor=1024,
              desc=f'      {fname}', leave=False) as pbar:
        with open(local_path, 'rb') as fh:
            wrapped = CallbackIOWrapper(pbar.update, fh, 'read')
            resp = requests.put(
                f'{bucket_url}/{fname}',
                data=wrapped,
                headers={
                    'Authorization':  f'Bearer {token}',
                    'Content-Type':   'application/octet-stream',
                    'Content-Length': str(size),
                },
            )
    if not resp.ok:
        raise requests.HTTPError(
            f'{resp.status_code} {resp.reason} — {resp.text}', response=resp)
    print(f'    ✓ {fname}', flush=True)


def publish_deposit(dep_id: str, token: str, dry_run: bool) -> None:
    url = f'{ZENODO_API}/deposit/depositions/{dep_id}/actions/publish'
    print(f'  → publishing deposit {dep_id}', flush=True)
    if dry_run:
        return
    resp = requests.post(url, headers={'Authorization': f'Bearer {token}'})
    resp.raise_for_status()
    print(f'  ✓ published', flush=True)


# ── Excel reader ─────────────────────────────────────────────────────────────

def load_upload_map(excel_path: str) -> dict[str, dict]:
    """
    Returns {upload_folder: {'doi': ..., 'dep_id': ..., 'files': [abbr, ...]}}
    for all redistributable studies that have a DOI assigned.
    """
    wb = openpyxl.load_workbook(excel_path)
    ws = wb['SupTab1']

    # Locate header row
    header_row_idx = None
    for i, row in enumerate(ws.iter_rows(values_only=True), 1):
        if row[0] == 'Trait_Number':
            header_row_idx = i
            break
    if header_row_idx is None:
        sys.exit('ERROR: could not find header row in SupTab1')

    uploads: dict[str, dict] = {}
    for row in ws.iter_rows(min_row=header_row_idx + 1, values_only=True):
        if not row[0]:
            continue
        abbr   = row[COL_ABBR]
        doi    = row[COL_DOI]
        folder = row[COL_UPLOAD]
        pub    = str(row[COL_PUB] or '')

        # Only redistributable studies with a DOI and upload folder
        if not doi or not folder or 'Free' not in pub:
            continue

        uploaded = bool(row[COL_UPLOADED])

        if folder not in uploads:
            uploads[folder] = {
                'doi':    doi,
                'dep_id': deposit_id_from_doi(doi),
                'files':  [],
            }
        uploads[folder]['files'].append({'abbr': abbr, 'uploaded': uploaded})

    return uploads


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--excel',     required=True,
                        help='Path to Tables.xlsx containing SupTab1')
    parser.add_argument('--per-study', required=True,
                        help='Path to the per_study directory with upload folders')
    parser.add_argument('--upload',    nargs='+', metavar='FOLDER',
                        help='Limit to these upload folder prefixes (e.g. upload17 upload18)')
    parser.add_argument('--publish',   action='store_true',
                        help='Publish each deposit after all its files are uploaded')
    parser.add_argument('--dry-run',   action='store_true',
                        help='Show what would be uploaded without transferring anything')
    args = parser.parse_args()

    token = os.environ.get('ZENODO_TOKEN')
    if not token:
        sys.exit('ERROR: ZENODO_TOKEN environment variable not set')

    per_study = Path(args.per_study)
    if not per_study.is_dir():
        sys.exit(f'ERROR: --per-study directory not found: {per_study}')

    print(f'Reading SupTab1 from {args.excel} …')
    upload_map = load_upload_map(args.excel)
    print(f'Found {len(upload_map)} upload records with DOIs.\n')

    # Filter to requested uploads
    if args.upload:
        prefixes = [u.strip('/') for u in args.upload]
        upload_map = {k: v for k, v in upload_map.items()
                      if any(k.startswith(p) for p in prefixes)}
        if not upload_map:
            sys.exit('ERROR: no matching upload folders found for --upload filter')

    total_uploads = 0
    total_skipped = 0

    for folder, info in sorted(upload_map.items()):
        dep_id     = info['dep_id']
        abbrs      = info['files']
        folder_dir = per_study / folder

        if not folder_dir.is_dir():
            print(f'[{folder}] WARNING: directory not found, skipping')
            continue

        n_studies  = len(abbrs)
        n_uploaded = sum(1 for f in abbrs if f['uploaded'])

        if n_uploaded == n_studies:
            print(f'[{folder}]  deposit={dep_id} — all {n_studies} studies already uploaded '
                  f'(column W), skipping')
            print()
            continue

        print(f'[{folder}]  deposit={dep_id}  '
              f'({n_studies} studies, {n_uploaded} already uploaded)')

        # Always fetch bucket state (read-only); only skip the PUT in dry-run.
        try:
            bucket_url = get_bucket_url(dep_id, token)
            existing   = list_bucket_files(bucket_url, token)
        except requests.HTTPError as e:
            print(f'  ERROR fetching deposit: {e}')
            continue

        # Collect files to upload: archives + manifest
        targets: list[Path] = []
        for file_info in abbrs:
            abbr = file_info['abbr']
            if file_info['uploaded']:
                print(f'  – {abbr}.tar.gz already uploaded (column W), skipping')
                total_skipped += 1
                continue
            p = folder_dir / f'{abbr}.tar.gz'
            if p.exists():
                targets.append(p)
            else:
                print(f'  WARNING: {p.name} not found')

        manifest = next(folder_dir.glob('manifest_upload*.tsv'), None)
        if manifest:
            targets.append(manifest)
        else:
            print(f'  WARNING: no manifest TSV found in {folder}')

        for fpath in targets:
            bucket_size = existing.get(fpath.name)
            if bucket_size is not None:
                local_size = fpath.stat().st_size
                if bucket_size == local_size:
                    print(f'    – {fpath.name} already in bucket, skipping')
                    total_skipped += 1
                    continue
                print(f'    ⚠ {fpath.name} partial in bucket '
                      f'({bucket_size / 1e9:.2f} / {local_size / 1e9:.2f} GB), re-uploading')
            upload_file(bucket_url, fpath, token, dry_run=args.dry_run)
            total_uploads += 1

        if args.publish and not args.dry_run:
            publish_deposit(dep_id, token, dry_run=False)
        elif args.publish and args.dry_run:
            print(f'  → [dry-run] would publish deposit {dep_id}')

        print()

    print(f'Done.  Uploaded: {total_uploads}  Skipped (already present): {total_skipped}')


if __name__ == '__main__':
    main()

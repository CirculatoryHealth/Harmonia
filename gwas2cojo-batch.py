#!/usr/bin/env python3
import argparse
import csv
import os
import shlex
import subprocess
import sys
from pathlib import Path
import time


# ----------------------------------------------------------------------
# Version info and copyright
# ----------------------------------------------------------------------
VERSION_NAME = "gwas2cojo_batch"
VERSION = "1.0.0"
VERSION_DATE = "2026-02-11"
COPYRIGHT = 'Copyright 1979-2026. Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
COPYRIGHT_TEXT = '''
The MIT License (MIT).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
'''

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"""
{VERSION_NAME} v{VERSION} ({VERSION_DATE})
Batch runner to execute gwas2cojo for multiple GWAS summary files listed in a CSV.""",
epilog=f"""
For a given list of GWAS files, run gwas2cojo.py for each file.
You must provide a CSV file with at least two columns:
- one column with the output name (e.g. CAD_aragam)
- one column with the GWAS file path

Optionally, you can specify the names of these columns via --output-col and --gwas-col.

Example:
    python gwas2cojo-batch.py -i gwas_list.csv --output-col output_name --gwas-col gwas_path -O ./cojo_results

+ {VERSION_NAME} v{VERSION}. {COPYRIGHT} +
{COPYRIGHT_TEXT}""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input-csv",
        required=True,
        help="Comma-separated file listing GWAS files and output names."
    )

    parser.add_argument(
        "--output-col",
        default="output_name",
        help="Name of the column containing the base output name (e.g. CAD_aragam)."
    )

    parser.add_argument(
        "--gwas-col",
        default="gwas_path",
        help="Name of the column containing the GWAS file path."
    )

    parser.add_argument(
        "--gwas-build-col",
        default="gwas_build",
        help="Optional column name in the CSV containing the GWAS build per row (e.g. hg19, hg38). If missing/empty, --gwas-build is used."
    )

    parser.add_argument(
        "--extra-args-col",
        default="extra_args",
        help="Optional column name in the CSV containing extra gwas2cojo CLI args per row (e.g. \"--gwas:build hg38 --some-flag\")."
    )

    parser.add_argument(
        "-O", "--output-dir",
        default="../../_cojo_2025",
        help="Base directory where -o and -r prefixes will be written."
    )

    parser.add_argument(
        "--gwas2cojo",
        default="/hpc/local/Rocky8/dhl_ec/software/gwas2cojo/gwas2cojo.py",
        help="Path to gwas2cojo.py."
    )

    parser.add_argument(
        "--gen",
        default="../../../references/1000G/Phase3/1kGp3v5b.ref.allfreq.noCN_noINS_noSS_noESV_noMultiAllelic.sumstats.txt.gz",
        help="Path to the reference --gen file."
    )

    parser.add_argument(
        "--gwas-build",
        default="hg19",
        help="Value for --gwas:build."
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands instead of executing them."
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output."
    )

    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"{VERSION_NAME} {VERSION} ({VERSION_DATE})"
    )

    return parser.parse_args()

# ----------------------------------------------------------------------
# Run command
# ----------------------------------------------------------------------
def run_command(cmd, dry_run: bool = False, verbose: bool = False) -> int:
    """Run a shell command via subprocess.run."""
    cmd_str = shlex.join(cmd)
    if dry_run or verbose:
        print(f"[CMD] {cmd_str}", flush=True)

    if dry_run:
        return 0

    result = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
    return result.returncode

# ----------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------
def main():
    
    script_start = time.time()

    args = parse_args()

    gwas2cojo_path = Path(args.gwas2cojo)
    if not gwas2cojo_path.is_file():
        print(f"ERROR: gwas2cojo.py not found: {gwas2cojo_path}", file=sys.stderr)
        sys.exit(1)

    gen_path = Path(args.gen)
    if not gen_path.is_file():
        print(f"ERROR: --gen reference file not found: {gen_path}", file=sys.stderr)
        sys.exit(1)

    input_csv = Path(args.input_csv)
    if not input_csv.is_file():
        print(f"ERROR: CSV file not found: {input_csv}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_total = 0
    n_run = 0
    n_skipped = 0
    n_failed = 0

    
    with input_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        missing_cols = [c for c in (args.output_col, args.gwas_col) if c not in (reader.fieldnames or [])]
        if missing_cols:
            print(
                f"ERROR: Missing expected columns in CSV: {', '.join(missing_cols)}\n"
                f"       Available columns: {', '.join(reader.fieldnames or [])}",
                file=sys.stderr,
            )
            sys.exit(1)

        for i, row in enumerate(reader, start=1):
            n_total += 1
            output_name = (row[args.output_col] or "").strip()
            gwas_path = (row[args.gwas_col] or "").strip()

            row_build = ""
            if args.gwas_build_col in row and row.get(args.gwas_build_col) is not None:
                row_build = str(row.get(args.gwas_build_col)).strip()
            gwas_build = row_build if row_build else args.gwas_build

            extra_args_raw = ""
            if args.extra_args_col in row and row.get(args.extra_args_col) is not None:
                extra_args_raw = str(row.get(args.extra_args_col)).strip()
            try:
                extra_args = shlex.split(extra_args_raw) if extra_args_raw else []
            except ValueError as e:
                print(f"[ERROR] Row {i}: could not parse {args.extra_args_col}='{extra_args_raw}': {e}", file=sys.stderr)
                n_failed += 1
                continue

            if not output_name:
                print(f"[WARN] Row {i}: empty output name, skipping.", file=sys.stderr)
                n_skipped += 1
                continue
            if not gwas_path:
                print(f"[WARN] Row {i}: empty GWAS path, skipping.", file=sys.stderr)
                n_skipped += 1
                continue

            gwas_path = os.path.expanduser(gwas_path)
            if not os.path.isabs(gwas_path):
                gwas_path = str((input_csv.parent / gwas_path).resolve())
            else:
                gwas_path = os.path.abspath(gwas_path)

            if not os.path.isfile(gwas_path):
                print(
                    f"[WARN] Row {i}: GWAS file does not exist: {gwas_path}. Skipping.",
                    file=sys.stderr,
                )
                n_skipped += 1
                continue

            # Build output prefix (same for -o and -r)
            out_prefix = output_dir / output_name

            cmd = [
                sys.executable,
                args.gwas2cojo,
                "-o", str(out_prefix),
                "-r", str(out_prefix),
                "--gen", args.gen,
                "--gwas", gwas_path,
                "--output-pos",
                "--gwas:build", gwas_build,
                "--gen:ident", "ID",
                "--gen:chr", "CHROM",
                "--gen:bp", "POS",
                "--gen:effect", "ALT",
                "--gen:other", "REF",
                "--gen:eaf", "EUR_AF",
            ]
            if extra_args:
                cmd.extend(extra_args)

            if args.verbose:
                print(f"\n[INFO] Row {i}: running gwas2cojo for output '{output_name}'")
                print(f"[INFO]   GWAS: {gwas_path}")
                print(f"[INFO]   Prefix: {out_prefix}")
                print(f"[INFO]   GWAS build: {gwas_build}")
                if extra_args:
                    print(f"[INFO]   Extra args: {shlex.join(extra_args)}")

            start_row = time.time()
            rc = run_command(cmd, dry_run=args.dry_run, verbose=args.verbose)
            elapsed_row = time.time() - start_row
            if args.verbose:
                print(f"[INFO]   Finished row {i} in {elapsed_row:.1f} seconds")

            if rc != 0:
                print(
                    f"[ERROR] gwas2cojo failed with exit code {rc} for row {i} (output '{output_name}').",
                    file=sys.stderr,
                )
                n_failed += 1
                continue

            n_run += 1

    if args.verbose:
        print(f"\n[INFO] All rows processed. total={n_total}, "
              f"run={n_run}, skipped={n_skipped}, failed={n_failed}")
        total_elapsed = time.time() - script_start
        print(f"[INFO] Total runtime: {total_elapsed:.1f} seconds")

# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
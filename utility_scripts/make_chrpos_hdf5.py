#!/usr/bin/env python3
"""
make_chrpos_hdf5.py — One-time conversion of dbSNP VCF to HDF5 for rsID→CHR:POS lookup.

Wraps gwaslab's gl.process_vcf_to_hfd5() to generate per-chromosome HDF5 files
(one per chromosome, 10 modulo groups each) used by gwas_process.py --add-chrpos.

This step needs to be run once per build. The HDF5 files are written into the
same REF_DIR as the VCF and are automatically discovered by gwas_process.py.

Usage
-----
    python make_chrpos_hdf5.py [--ref-dir DIR] [--build hg19|hg38|all]
                               [--threads N] [--complevel 0-9] [--overwrite]

    # Use REF_DIR from gwas2cojo.conf (default):
    python make_chrpos_hdf5.py --build hg19

    # Override reference directory:
    python make_chrpos_hdf5.py --ref-dir /hpc/data/references/gwaslab --build hg19

    # Regenerate (overwrite existing HDF5 files):
    python make_chrpos_hdf5.py --build hg19 --overwrite

Requirements
------------
    - gwaslab >= 3.4.42 (process_vcf_to_hfd5 available)
    - bcftools in PATH (used internally by gwaslab)
    - dbSNP VCF already downloaded (run gwas_process.download_refs.py first)

Output
------
    {REF_DIR}/GCF_*.chr{N}.rsID_CHR_POS_mod10.h5   (one per chromosome)

Runtime
-------
    Approximately 1-3 hours per build on 6 threads for a ~24 GB dbSNP VCF.
"""

# ── Version ───────────────────────────────────────────────────────────────────
VERSION_NAME = "make_chrpos_hdf5"
VERSION      = "1.0.0"
VERSION_DATE = "2026-04-03"

import argparse
import os
import re
import sys

# ── Read REF_DIR from gwas2cojo.conf ─────────────────────────────────────────
_FALLBACK_REF_DIR = "/path/to/references/gwaslab/"

def _read_conf_ref_dir() -> str:
    conf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "..", "gwas2cojo.conf")
    conf_path = os.path.normpath(conf_path)
    if not os.path.isfile(conf_path):
        return _FALLBACK_REF_DIR
    with open(conf_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or "=" not in line:
                continue
            m = re.match(r'^REF_DIR=["\']?([^"\'#\n]+?)["\']?\s*(?:#.*)?$', line)
            if m:
                return m.group(1).strip()
    return _FALLBACK_REF_DIR


DEFAULT_REF_DIR = _read_conf_ref_dir()

# ── dbSNP VCF gwaslab keys per build ─────────────────────────────────────────
# gwas_process.download_refs.py downloads these; we prefer v157 (newer), fall
# back to v151 if v157 is not present.
_DBSNP_KEYS = {
    "hg19": ["dbsnp_v157_hg19", "dbsnp_v151_hg19"],
    "hg38": ["dbsnp_v157_hg38", "dbsnp_v151_hg38"],
}


def _find_vcf(ref_dir: str, build: str) -> str:
    """Return path to the best available dbSNP VCF for the given build."""
    try:
        import gwaslab as gl
        from gwaslab.bd.bd_download import set_default_directory, get_path
        set_default_directory(ref_dir)
    except ImportError:
        print("ERROR: gwaslab is not installed.", file=sys.stderr)
        sys.exit(1)

    for key in _DBSNP_KEYS[build]:
        try:
            path = get_path(key, verbose=False)
            if path and os.path.isfile(path):
                print(f"  Found dbSNP VCF ({key}): {path}")
                return path
        except Exception:
            continue

    print(f"ERROR: no dbSNP VCF found for {build} in {ref_dir}", file=sys.stderr)
    print(f"       Run: python gwas_process.download_refs.py --build {build[2:]}", file=sys.stderr)
    sys.exit(1)


def make_hdf5(ref_dir: str, build: str, threads: int, complevel: int,
              overwrite: bool) -> None:
    """Generate HDF5 files for one build."""
    import gwaslab as gl

    vcf_path = _find_vcf(ref_dir, build)

    print(f"\n{'=' * 60}")
    print(f"Build        : {build}")
    print(f"VCF          : {vcf_path}")
    print(f"Output dir   : {ref_dir}")
    print(f"Threads      : {threads}")
    print(f"Compression  : {complevel}")
    print(f"Overwrite    : {overwrite}")
    print(f"{'=' * 60}\n")

    gl.process_vcf_to_hfd5(
        vcf       = vcf_path,
        directory = ref_dir,
        complevel = complevel,
        threads   = threads,
        overwrite = overwrite,
    )
    print(f"\nHDF5 files written to: {ref_dir}")


def main() -> None:
    p = argparse.ArgumentParser(
        prog="make_chrpos_hdf5.py",
        description=(
            f"make_chrpos_hdf5.py  v{VERSION}  ({VERSION_DATE})\n"
            "One-time conversion of dbSNP VCF → per-chromosome HDF5 files\n"
            "for fast rsID→CHR:POS lookup in gwas_process.py --add-chrpos."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--ref-dir",   default=DEFAULT_REF_DIR, metavar="DIR",
                   help=f"gwaslab reference directory (default: {DEFAULT_REF_DIR}).")
    p.add_argument("--build",     default="hg19",
                   choices=["hg19", "hg38", "all"],
                   help="Genome build to process (default: hg19). "
                        "Use 'all' to process both hg19 and hg38.")
    p.add_argument("--threads",   default=6, type=int, metavar="N",
                   help="Parallel threads for processing (default: 6).")
    p.add_argument("--complevel", default=3, type=int, metavar="0-9",
                   help="HDF5 compression level 0-9 (default: 3).")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing HDF5 files (default: skip existing).")
    args = p.parse_args()

    if not os.path.isdir(args.ref_dir):
        print(f"ERROR: ref-dir does not exist: {args.ref_dir}", file=sys.stderr)
        sys.exit(1)

    builds = ["hg19", "hg38"] if args.build == "all" else [args.build]
    for build in builds:
        make_hdf5(
            ref_dir   = args.ref_dir,
            build     = build,
            threads   = args.threads,
            complevel = args.complevel,
            overwrite = args.overwrite,
        )

    print("\nDone. HDF5 files are ready for use with gwas_process.py --add-chrpos.")


if __name__ == "__main__":
    main()

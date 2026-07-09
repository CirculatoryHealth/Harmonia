#!/usr/bin/env python3
"""
make_chrpos_hdf5.py — One-time conversion of dbSNP VCF to HDF5 for rsID→CHR:POS lookup.

Wraps gwaslab's gl.process_vcf_to_hfd5() to generate per-chromosome HDF5 files
(one per chromosome, 10 modulo groups each) used by harmonia.py --add-chrpos.

This step needs to be run once per build. The HDF5 files are written into the
same REF_DIR as the VCF and are automatically discovered by harmonia.py.

Usage
-----
    python make_chrpos_hdf5.py [--ref-dir DIR] [--build hg19|hg38|all]
                               [--threads N] [--complevel 0-9] [--overwrite]
                               [--vcf /path/to/dbsnp.vcf.gz]

    # Use REF_DIR from gwas2cojo.conf (default):
    python make_chrpos_hdf5.py --build hg19

    # Override reference directory:
    python make_chrpos_hdf5.py --ref-dir /hpc/data/references/gwaslab --build hg19

    # Specify the VCF directly (when auto-detection fails, e.g. GCF-named files):
    python make_chrpos_hdf5.py --build hg19 \
        --vcf /hpc/data/references/gwaslab/GCF_000001405.25.gz

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
VERSION      = "1.1.0"
VERSION_DATE = "2026-04-04"

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

# ── GRCh38 RefSeq accession → chromosome number mapping ──────────────────────
# dbSNP b157 VCFs for GRCh38 have no ##contig header lines. gwaslab's
# auto-detection therefore cannot distinguish hg38 from hg19 RefSeq IDs and
# falls back to the hg19 mapping (NC_000001.10 → 1, etc.).  Passing this dict
# explicitly tells bcftools to look for the hg38 accessions (NC_000001.11 etc.)
# instead, which are the actual CHROM values in the GRCh38 data records.
_REFSEQ_HG38 = {
    "NC_000001.11": "1",  "NC_000002.12": "2",  "NC_000003.12": "3",
    "NC_000004.12": "4",  "NC_000005.10": "5",  "NC_000006.12": "6",
    "NC_000007.14": "7",  "NC_000008.11": "8",  "NC_000009.12": "9",
    "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9":  "14", "NC_000015.10": "15",
    "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
    "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9":  "21",
    "NC_000022.11": "22", "NC_000023.11": "23", "NC_000024.10": "24",
    "NC_012920.1":  "25",
}

# ── dbSNP VCF gwaslab keys per build ─────────────────────────────────────────
# gwas_process.download_refs.py downloads these; we prefer v157 (newer), fall
# back to v151 if v157 is not present.
_DBSNP_KEYS = {
    "hg19": ["dbsnp_v157_hg19", "dbsnp_v151_hg19"],
    "hg38": ["dbsnp_v157_hg38", "dbsnp_v151_hg38"],
}

# ── GCF accession prefixes per build ─────────────────────────────────────────
# NCBI names dbSNP VCF files after the RefSeq assembly accession.
# GCF_000001405.25 = GRCh37.p13 (hg19)
# GCF_000001405.40 = GRCh38.p14 (hg38, latest)
# GCF_000001405.39 = GRCh38.p13 (hg38, older)
_DBSNP_GCF_PREFIXES = {
    "hg19": ["GCF_000001405.25"],
    "hg38": ["GCF_000001405.40", "GCF_000001405.39"],
}


def _find_vcf(ref_dir: str, build: str, vcf_override: str = "") -> str:
    """Return path to the best available dbSNP VCF for the given build.

    Detection order:
      1. ``--vcf`` direct override (skips all auto-detection).
      2. gwaslab key lookup (dbsnp_v157_*/dbsnp_v151_*) via get_path().
      3. GCF accession glob in ref_dir (GCF_000001405.25*.gz for hg19, etc.).
      4. Classic NCBI name ``00-All.vcf.gz`` in ref_dir (build-ambiguous;
         logged as a warning).
    """
    # ── 1. Direct path override ───────────────────────────────────────────────
    if vcf_override:
        if os.path.isfile(vcf_override):
            print(f"  Using --vcf override: {vcf_override}")
            return vcf_override
        print(f"ERROR: --vcf file not found: {vcf_override}", file=sys.stderr)
        sys.exit(1)

    # ── 2. gwaslab key lookup ─────────────────────────────────────────────────
    try:
        from gwaslab.bd.bd_download import set_default_directory, get_path
        set_default_directory(ref_dir)
        for key in _DBSNP_KEYS[build]:
            try:
                path = get_path(key, verbose=False)
                if path and os.path.isfile(path):
                    print(f"  Found dbSNP VCF (gwaslab key '{key}'): {path}")
                    return path
            except Exception:
                continue
    except ImportError:
        print("ERROR: gwaslab is not installed.", file=sys.stderr)
        sys.exit(1)

    # ── 3. GCF accession pattern fallback ────────────────────────────────────
    # Files in ref_dir downloaded by NCBI are named after the assembly accession
    # (e.g. GCF_000001405.25.gz) rather than gwaslab's internal key names.
    for prefix in _DBSNP_GCF_PREFIXES.get(build, []):
        try:
            candidates = [
                f for f in os.listdir(ref_dir)
                if f.startswith(prefix) and f.endswith(".gz")
                and not f.endswith(".gz.tbi")
            ]
        except OSError:
            candidates = []
        for fname in sorted(candidates):          # deterministic: pick shortest/oldest
            path = os.path.join(ref_dir, fname)
            print(f"  Found dbSNP VCF (GCF pattern '{prefix}'): {path}")
            return path

    # ── 4. Classic 00-All.vcf.gz fallback ────────────────────────────────────
    fallback = os.path.join(ref_dir, "00-All.vcf.gz")
    if os.path.isfile(fallback):
        print(f"  WARNING: using 00-All.vcf.gz — verify this file is {build}.",
              file=sys.stderr)
        print(f"  Found dbSNP VCF (00-All fallback): {fallback}")
        return fallback

    print(f"ERROR: no dbSNP VCF found for {build} in {ref_dir}", file=sys.stderr)
    print(f"       Tried gwaslab keys: {_DBSNP_KEYS[build]}", file=sys.stderr)
    print(f"       Tried GCF prefixes: {_DBSNP_GCF_PREFIXES.get(build, [])}", file=sys.stderr)
    print(f"       Tried fallback    : {fallback}", file=sys.stderr)
    print(f"       Use --vcf /path/to/dbsnp_{build}.vcf.gz to specify directly.",
          file=sys.stderr)
    sys.exit(1)


def make_hdf5(ref_dir: str, build: str, threads: int, complevel: int,
              overwrite: bool, vcf_override: str = "") -> None:
    """Generate HDF5 files for one build."""
    import gwaslab as gl

    vcf_path = _find_vcf(ref_dir, build, vcf_override=vcf_override)

    print(f"\n{'=' * 60}")
    print(f"Build        : {build}")
    print(f"VCF          : {vcf_path}")
    print(f"Output dir   : {ref_dir}")
    print(f"Threads      : {threads}")
    print(f"Compression  : {complevel}")
    print(f"Overwrite    : {overwrite}")
    print(f"{'=' * 60}\n")

    # For hg38, pass an explicit chr_dict so gwaslab uses the correct GRCh38
    # RefSeq accession IDs. Without it, gwaslab auto-detection falls back to
    # hg19 IDs (NC_000001.10 etc.) when the VCF has no ##contig header lines,
    # causing bcftools to find 0 rows for every autosome.
    chr_dict = _REFSEQ_HG38 if build == "hg38" else None

    gl.process_vcf_to_hfd5(
        vcf       = vcf_path,
        directory = ref_dir,
        complevel = complevel,
        threads   = threads,
        overwrite = overwrite,
        chr_dict  = chr_dict,
    )
    print(f"\nHDF5 files written to: {ref_dir}")


def main() -> None:
    p = argparse.ArgumentParser(
        prog="make_chrpos_hdf5.py",
        description=(
            f"make_chrpos_hdf5.py  v{VERSION}  ({VERSION_DATE})\n"
            "One-time conversion of dbSNP VCF → per-chromosome HDF5 files\n"
            "for fast rsID→CHR:POS lookup in harmonia.py --add-chrpos."
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
    p.add_argument("--vcf",      default="", metavar="PATH",
                   help="Direct path to the dbSNP VCF (.vcf.gz).  Bypasses automatic "
                        "detection (gwaslab key lookup + GCF pattern scan).  Useful "
                        "when the file is named differently from gwaslab's conventions "
                        "(e.g. GCF_000001405.25.gz or 00-All.vcf.gz).  "
                        "Not compatible with --build all.")
    args = p.parse_args()

    if args.vcf and args.build == "all":
        print("ERROR: --vcf cannot be combined with --build all "
              "(--vcf specifies a single file for a single build).", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.ref_dir):
        print(f"ERROR: ref-dir does not exist: {args.ref_dir}", file=sys.stderr)
        sys.exit(1)

    builds = ["hg19", "hg38"] if args.build == "all" else [args.build]
    for build in builds:
        make_hdf5(
            ref_dir      = args.ref_dir,
            build        = build,
            threads      = args.threads,
            complevel    = args.complevel,
            overwrite    = args.overwrite,
            vcf_override = args.vcf,
        )

    print("\nDone. HDF5 files are ready for use with harmonia.py --add-chrpos.")


if __name__ == "__main__":
    main()

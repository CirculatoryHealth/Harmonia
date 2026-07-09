#!/usr/bin/env python3
"""
gwas_process.download_refs.py — Download Harmonia/GWASLab reference files.

Usage
-----
    python gwas_process.download_refs.py [--ref-dir DIR] [--build all|hg19|hg38]
                                    [--ancestry EUR|AFR|EAS|AMR|SAS|PAN|all]
                                    [--no-x]

The reference directory defaults to REF_DIR from gwas2cojo.conf (looked up
next to this script).  Pass --ref-dir to override.

--build controls which genome build(s) to download VCFs and build-specific
files for (default: all).

--ancestry controls which 1KG population VCF(s) to download (default: EUR).
Use 'all' to download every available ancestry.

Chromosome X VCFs (1kg_*_x_hg19/hg38) and chrX SNPID→rsID tables
(1kg_dbsnp151_hg19_x / hg38_x) are downloaded by default alongside the
autosomal files.  Pass --no-x to skip all chrX reference files.

Non-ancestry-specific files (dbSNP VCFs, FASTA, recombination maps, GTFs,
HapMap3 EAF, SNPID→rsID tables) are always included for the requested build(s).

Note: recombination maps and GTF files are downloaded automatically by gwaslab
at runtime when needed; pre-downloading them here is optional but useful on HPC
nodes without outbound internet during jobs.

Reference
---------
https://cloufield.github.io/gwaslab/Download/
"""

# ============================================================
VERSION_NAME = "harmonia.download_refs"
VERSION      = "1.2.1"
VERSION_DATE = "2026-03-27"
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

# ── Loading packages ──────────────────────────────────────────────────────────
import argparse
import os
import re
import gwaslab as gl
from gwaslab.bd.bd_download import (
    set_default_directory,
    get_default_directory,
    download_ref,
    check_downloaded_ref,
    get_path,
)

_FALLBACK_REF_DIR = "/path/to/references/gwaslab/"

ALL_ANCESTRIES = ["EUR", "PAN", "AFR", "EAS", "AMR", "SAS"]
ALL_BUILDS     = ["hg19", "hg38"]

# ───────────────────────────────────────────────────────────────────────────
# Reference file catalogue — organised by (ancestry, build).
# ancestry=None means the file is not population-specific.
# ───────────────────────────────────────────────────────────────────────────
# fmt: off
_ANCESTRY_VCFS = {
    # 1KG population VCFs — hg19 (1KGp3v5 low-coverage)
    ("EUR", "hg19"): "1kg_eur_hg19",
    ("PAN", "hg19"): "1kg_pan_hg19",
    ("AFR", "hg19"): "1kg_afr_hg19",
    ("EAS", "hg19"): "1kg_eas_hg19",
    ("AMR", "hg19"): "1kg_amr_hg19",
    ("SAS", "hg19"): "1kg_sas_hg19",
    # 1KG population VCFs — hg38 (1KGp 30x)
    ("EUR", "hg38"): "1kg_eur_hg38",
    ("PAN", "hg38"): "1kg_pan_hg38",
    ("AFR", "hg38"): "1kg_afr_hg38",
    ("EAS", "hg38"): "1kg_eas_hg38",
    ("AMR", "hg38"): "1kg_amr_hg38",
    ("SAS", "hg38"): "1kg_sas_hg38",
}

_ANCESTRY_X_VCFS = {
    # 1KG chrX VCFs — hg19 (1KGp3v5)
    ("EUR", "hg19"): "1kg_eur_x_hg19",
    ("PAN", "hg19"): "1kg_pan_x_hg19",
    ("AFR", "hg19"): "1kg_afr_x_hg19",
    ("EAS", "hg19"): "1kg_eas_x_hg19",
    ("AMR", "hg19"): "1kg_amr_x_hg19",
    ("SAS", "hg19"): "1kg_sas_x_hg19",
    # 1KG chrX VCFs — hg38 (1KGp 30x)
    ("EUR", "hg38"): "1kg_eur_x_hg38",
    ("PAN", "hg38"): "1kg_pan_x_hg38",
    ("AFR", "hg38"): "1kg_afr_x_hg38",
    ("EAS", "hg38"): "1kg_eas_x_hg38",
    ("AMR", "hg38"): "1kg_amr_x_hg38",
    ("SAS", "hg38"): "1kg_sas_x_hg38",
}

_BUILD_FILES = {
    "hg19": [
        "1kg_hm3_hg19_eaf",        # HapMap3 pan-ancestry EAF
        "1kg_dbsnp151_hg19_auto",  # SNPID→rsID conversion table (autosomes)
        "1kg_dbsnp151_hg19_x",     # SNPID→rsID conversion table (chrX)
        "dbsnp_v151_hg19",         # dbSNP v151 VCF  (!large)
        "dbsnp_v157_hg19",         # dbSNP v157 VCF  (!large)
        "ucsc_genome_hg19",        # reference FASTA  (!large)
        "recombination_hg19",      # recombination map
        "ensembl_hg19_gtf",        # Ensembl gene annotation GTF
        "refseq_hg19_gtf",         # RefSeq gene annotation GTF
    ],
    "hg38": [
        "1kg_hm3_hg38_eaf",
        "1kg_dbsnp151_hg38_auto",  # SNPID→rsID conversion table (autosomes)
        "1kg_dbsnp151_hg38_x",     # SNPID→rsID conversion table (chrX)
        "dbsnp_v151_hg38",
        "dbsnp_v157_hg38",
        "ucsc_genome_hg38",
        "recombination_hg38",
        "ensembl_hg38_gtf",
        "refseq_hg38_gtf",
    ],
}
# fmt: on

# ── Read configuration ───────────────────────────────────────────────────────────
def _read_conf_ref_dir() -> str:
    """
    Read REF_DIR from gwas2cojo.conf located next to this script.
    Returns the fallback placeholder string if the conf is absent or
    REF_DIR is not set.
    """
    conf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gwas2cojo.conf")
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

# ── Download list function ───────────────────────────────────────────────────────────
def build_download_list(builds: list, ancestries: list, include_x: bool = True) -> list:
    """Return the ordered list of gwaslab keywords to download."""
    to_download = []
    # 1KG ancestry VCFs — autosomes first, then chrX
    for ancestry in ancestries:
        for build in builds:
            key = (ancestry.upper(), build)
            if key in _ANCESTRY_VCFS:
                to_download.append(_ANCESTRY_VCFS[key])
    if include_x:
        for ancestry in ancestries:
            for build in builds:
                key = (ancestry.upper(), build)
                if key in _ANCESTRY_X_VCFS:
                    to_download.append(_ANCESTRY_X_VCFS[key])
    # Build-specific but ancestry-agnostic files
    # (chrX SNPID→rsID tables are included in _BUILD_FILES when include_x is True;
    #  skip them here if --no-x was requested)
    for build in builds:
        files = _BUILD_FILES[build]
        if not include_x:
            files = [f for f in files if "_x" not in f]
        to_download.extend(files)
    return to_download

# ── Main function ───────────────────────────────────────────────────────────
def main() -> None:
    # ── Argument parsing ───────────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        description=(
            f"gwas_process.download_refs  v{VERSION}  ({VERSION_DATE})\n"
            "Download gwaslab reference files for the gwas2cojo pipeline.\n"
            "Target directory is read from gwas2cojo.conf (REF_DIR) by default."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # EUR VCFs for both builds (default)\n"
            "  python gwas_process.download_refs.py\n\n"
            "  # EUR VCFs, hg38 only\n"
            "  python gwas_process.download_refs.py --build hg38\n\n"
            "  # AFR VCFs for both builds\n"
            "  python gwas_process.download_refs.py --ancestry AFR\n\n"
            "  # All ancestries, hg19 only\n"
            "  python gwas_process.download_refs.py --build hg19 --ancestry all\n\n"
            "  # Everything, custom directory\n"
            "  python gwas_process.download_refs.py --build all --ancestry all \\\n"
            "      --ref-dir /path/to/references/gwaslab/\n\n"
            "  # Autosomes only, skip chrX VCFs and chrX rsID tables\n"
            "  python gwas_process.download_refs.py --no-x\n"
        ),
    )
    parser.add_argument(
        "--ref-dir", default=DEFAULT_REF_DIR,
        help=f"Reference directory (default: read from gwas2cojo.conf, "
             f"currently: {DEFAULT_REF_DIR}).",
    )
    parser.add_argument(
        "--build", default="all", choices=["all", "hg19", "hg38"],
        help="Genome build(s) to download files for (default: all).",
    )
    parser.add_argument(
        "--ancestry", default="EUR", choices=["PAN", "AFR", "EAS", "AMR", "SAS", "EUR"],
        help=(
            "Ancestry/population VCF(s) to download. "
            f"Choices: {', '.join(ALL_ANCESTRIES)}, all  (default: EUR)."
        ),
    )
    parser.add_argument(
        "--no-x", dest="include_x", action="store_false", default=True,
        help=(
            "Skip chromosome X VCFs and chrX SNPID→rsID tables "
            "(by default chrX references are downloaded alongside autosomes)."
        ),
    )
    args = parser.parse_args()

    # Resolve ancestry list
    if args.ancestry.lower() == "all":
        ancestries = ALL_ANCESTRIES
    else:
        ancestry = args.ancestry.upper()
        if ancestry not in ALL_ANCESTRIES:
            parser.error(
                f"Unknown ancestry '{args.ancestry}'. "
                f"Choose from: {', '.join(ALL_ANCESTRIES)}, all"
            )
        ancestries = [ancestry]

    # Resolve build list
    builds = ALL_BUILDS if args.build == "all" else [args.build]

    # Print configuration
    if args.ref_dir == _FALLBACK_REF_DIR:
        print("WARNING: gwas2cojo.conf not found or REF_DIR not set.")
        print(f"         Using fallback: {_FALLBACK_REF_DIR}")
        print("         Copy gwas2cojo.conf.example → gwas2cojo.conf and set REF_DIR.")
    else:
        print(f"Reference directory : {args.ref_dir}  (from gwas2cojo.conf)")
    print(f"Build(s)            : {', '.join(builds)}")
    print(f"Ancestry/population : {', '.join(ancestries)}")
    print(f"Include chrX        : {'yes' if args.include_x else 'no (--no-x)'}")

    to_download = build_download_list(builds, ancestries, include_x=args.include_x)
    print(f"Files to download   : {len(to_download)}")

    set_default_directory(args.ref_dir)

    for keyword in to_download:
        print(f"\n{'─' * 60}")
        print(f"Downloading: {keyword}")
        download_ref(keyword)
        path = get_path(keyword, verbose=False)
        print(f"  → {path}")

    print(f"\n{'═' * 60}")
    print("Done. Downloaded references:")
    check_downloaded_ref()

# ───────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()

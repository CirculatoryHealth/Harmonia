#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# GWASLab Cohort Processing Pipeline
# ============================================================
# Standalone GWASLab processing pipeline for GWAS summary statistics.
#
# Default pipeline (all steps ON unless disabled with --no-* flags):
#
#   1.  Load data
#   2.  Verify / auto-detect genome build  [always]
#   3.  EAF check & optional lookup from 1KG reference VCF  [always]
#   4.  Column correction (P, SE, N derivation)  [always]
#   5.  Column standardisation  [always]
#   6.  Create GWASLab Sumstats object  [always]
#   7.  Plot raw histograms  [--no-figures to disable]
#   8.  basic_check  [always]
#   9.  remove_dup  [always]
#  10.  Liftover to hg38 if build is hg18 or hg19  [--liftover to enable]
#        hg19→hg38: built-in chain  |  hg18→hg38: requires hg18ToHg38.over.chain.gz in --ref
#  11.  check_ref + flip_allele_stats  [always, requires FASTA]
#  12.  fix_id  [always]
#  13.  infer_strand + flip_allele_stats  [always]
#  14.  assign_rsid (dbSNP)  [--no-dbsnp to disable]
#  15.  check_af2 [always]  (sweep variant — bcftools one-pass)
#  16.  Save raw output (pickle + parquet + TSV.GZ)  [always; use --no-pickle to skip pkl]
#  17.  Manhattan + QQ plots (unfiltered)  [--no-figures to disable]
#  18.  QC filter  [--no-qc to disable]
#  19.  Save QC output  [when QC enabled]
#  20.  Manhattan + QQ plots (QC-filtered)  [--no-figures to disable]
#  21.  Lead SNPs  [--no-leads to disable]
#  22.  COJO output  [--cojo to enable]
#
# COJO format:  SNP  A1  A2  freq  b  se  p  n
#   --cojo-id   chrpos (default) | rsid
#   --cojo-pos  include CHR and BP columns after SNP column
#
# Usage examples:
#
#   # Minimal — run full default pipeline
#   python3 gwaslab_process.py \
#       --gwas    MyStudy \
#       --input   MyStudy.parsed.txt.gz \
#       --dir     /data/results/MyStudy \
#       --ref     /data/references/ \
#       --pop     EUR \
#       --build   19
#
#   # Disable liftover and QC; produce COJO with rsID SNP column + CHR/BP cols
#   python3 gwaslab_process.py \
#       --gwas    MyStudy \
#       --input   MyStudy.parsed.txt.gz \
#       --dir     /data/results/MyStudy \
#       --ref     /data/references/ \
#       --pop     EUR \
#       --build   38 \
#       --no-liftover --no-qc \
#       --cojo --cojo-id rsid --cojo-pos


# ============================================================
VERSION_NAME = "gwaslab_process"
VERSION      = "1.4.1"
VERSION_DATE = "2026-03-18"
COPYRIGHT = 'Copyright 1979-2026. Emma J.A. Smulders; Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | https://vanderlaanand.science.'
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
# General-purpose packages
import sys
import os
import json
import argparse
import gzip
import logging
import shutil
import gc

# Visualisation and data handling
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for script use
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.stats import norm

# GWASLab
import gwaslab as gl

# ── Logging setup ─────────────────────────────────────────────────────────────

def setup_logging(log_path: str) -> logging.Logger:
    """Initialise root logger writing to both console and *log_path*."""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    # File handler
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="gwaslab_process.py",
        description=(
            f"GWASLab Cohort Processing Pipeline  v{VERSION}  ({VERSION_DATE})\n"
            "Standardise, harmonise, and QC GWAS summary statistics."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python3 gwaslab_process.py \\\n"
            "      --gwas MyStudy \\\n"
            "      --input MyStudy.parsed.txt.gz \\\n"
            "      --directory /data/results/MyStudy \\\n"
            "      --ref /data/references/ \\\n"
            "      --output /data/results/MyStudy/GWASCatalog \\\n"
            "      --population EUR \\\n"
            "      --build 19 \\\n"
            "      --dbsnp --qc\n"
        ),
    )

    # ── Study / paths ───────────────────────────────────────────────────────────
    req = p.add_argument_group("required arguments")
    req.add_argument("--gwas",  required=True, metavar="NAME",
                     help="Study / phenotype name; used in all output filenames.")
    req.add_argument("--input", required=True, metavar="FILE",
                     help="Input summary-statistics file (TSV or TSV.GZ). "
                          "Relative paths are resolved against --dir.")
    req.add_argument("--ref",   required=True, metavar="DIR",
                     help="Reference-files directory (FASTA, 1KG VCFs, dbSNP VCFs).")

    path = p.add_argument_group("path arguments")
    path.add_argument("--directory",    default=".", metavar="DIR",
                      dest="directory",
                      help="Directory containing the input file (default: '.').")
    path.add_argument("--output", default=None, metavar="DIR",
                      help="Output directory. If omitted, constructed at runtime as "
                           "<--directory>/<--gwas>/GWASCatalog (e.g. --directory /data and "
                           "--gwas MyStudy → /data/MyStudy/GWASCatalog). "
                           "The default cannot be shown here because it depends on "
                           "--directory and --gwas which are only known at runtime.")

    # ── Population / build ───────────────────────────────────────────────────────
    pop = p.add_argument_group("population / build")
    pop.add_argument("--population", default="EUR",
                     dest="population",
                     choices=["EUR", "AFR", "EAS", "AMR", "SAS", "PAN"],
                     help="Ancestry/population code used for 1KG reference VCF "
                          "(default: EUR).")
    pop.add_argument("--build", default="19",
                     choices=["18", "19", "38",
                              "hg18", "hg19", "hg38",
                              "GRCh36", "GRCh37", "GRCh38",
                              "b36", "b37", "b38"],
                     help="Input genome build (default: 19 = GRCh37). "
                          "Build 18/hg18/b36 triggers hg18→hg38 liftover using "
                          "hg18ToHg38.over.chain.gz from --ref.")

    # ── Pipeline toggles ─────────────────────────────────────────────────────────
    tog = p.add_argument_group("pipeline step toggles")
    tog.add_argument("--liftover",     action="store_true",
                   help="Lift over coordinates to the other build (hg19→hg38 or hg38→hg19).")
    tog.add_argument("--dbsnp",        action="store_true",
                   help="Assign rsIDs from a dbSNP VCF.")
    tog.add_argument("--qc",           action="store_true",
                   help="Apply QC filters and save a filtered output.")
    tog.add_argument("--only-qc",      action="store_true",
                   help="Skip loading; reload from an existing pickle and run QC only.")
    tog.add_argument("--fill-eaf",     action="store_true",
                   help="Look up missing EAF values from the (1KG) reference VCF "
                        "(slow for large datasets; skipped by default).")
    tog.add_argument("--no-pickle",    action="store_true",
                   help="Skip saving .pkl files (reduces peak memory and disk usage "
                        "on the save step; disables --only-qc for this run).")
    tog.add_argument("--stage",
                   choices=["all",
                            "preprocess",
                            "process-normalize",
                            "process-split",
                            "process-check-ref",
                            "process-infer-strand",
                            "process-assign-rsid",
                            "process-check-af",
                            "merge",
                            "qc",
                            "cojo"],
                   default="all", metavar="STAGE",
                   help=("Run a single pipeline stage using parquet/pickle checkpoints "
                         "for handoff, so each stage can be submitted as a separate "
                         "SLURM job with its own memory and time limits. "
                         "  preprocess          : load + standardise → .preprocess.parquet  (light)  "
                         "  process-normalize   : basic_check + remove_dup + liftover → .normalize.pkl  (medium)  "
                         "  process-split       : split normalize.pkl into per-chr parquets + manifest  (trivial)  "
                         "  process-check-ref   : check_ref + flip + fix_id  (medium; per-chr with --chrom)  "
                         "  process-infer-strand: infer_strand2 + flip  (high — 1KG sweep; per-chr with --chrom)  "
                         "  process-assign-rsid : assign_rsid  (extreme — dbSNP sweep; per-chr with --chrom)  "
                         "  process-check-af    : check_af2  (high — 1KG sweep; per-chr with --chrom)  "
                         "  merge               : concat per-chr parquets + QC + plots + leads + COJO  (medium)  "
                         "  qc                  : QC filter → .qc.pkl + .qc.parquet + .qc.tsv.gz  (medium; whole-genome only)  "
                         "  cojo                : write COJO file from existing pickle  (light; whole-genome only)  "
                         "  all                 : full pipeline end-to-end, no checkpoints (default). "
                         "Pass identical --gwas / --build / --liftover / --output / --dbsnp flags "
                         "to every stage so file stems and checkpoint paths match."))
    tog.add_argument("--chrom", type=int, default=None, metavar="N",
                   help=("Chromosome task ID 1–26 for per-chromosome processing. "
                         "23=X  24=Y  25=nonPAR  26=MT. "
                         "Set automatically by SLURM array jobs via $SLURM_ARRAY_TASK_ID. "
                         "When provided, process-check-ref / process-infer-strand / "
                         "process-assign-rsid / process-check-af each operate on a single "
                         "chromosome shard written by process-split. "
                         "If the chromosome is not present in the study the stage exits "
                         "gracefully (exit 0) so SLURM afterok dependencies are satisfied."))
    tog.add_argument("--figures",      action="store_true",
                   help="Generate diagnostic plots.")
    tog.add_argument("--leads",        action="store_true",
                   help="Extract genome-wide significant lead SNPs (p < 5e-8).")
    
    # ── COJO output ────────────────────────────────────────────────────────────
    cojo = p.add_argument_group("COJO output (gcta --cojo format)")
    cojo.add_argument("--cojo", action="store_true",
                      help="Write a COJO-format file "
                           "(columns: SNP A1 A2 freq b se p n).")
    cojo.add_argument("--cojo-id",
                      choices=["chrpos", "rsid"], default="chrpos",
                      help="SNPID format for the SNP column: "
                           "'chrpos' = CHR:BP (default), "
                           "'rsid'   = rsID (requires --dbsnp or pre-existing rsID).")
    cojo.add_argument("--cojo-pos", action="store_true",
                      help="Add CHR and BP columns immediately after the SNP column.")

    # ── Performance ───────────────────────────────────────────────────────────
    perf = p.add_argument_group("performance")
    perf.add_argument("--threads", type=int, default=4, metavar="N",
                      help="Threads for parallelised gwaslab steps (default: 4).")

    # ── Fixed sample sizes ───────────────────────────────────────────────────────
    nsz = p.add_argument_group(
        "fixed sample sizes",
        "Override or supply N columns when they are absent from the input file. "
        "Values from the paper can be provided here. If the column already exists "
        "in the data it is left unchanged (use --force-n to override).",
    )
    nsz.add_argument("--n",          type=int, default=None, metavar="INT",
                     help="Total sample size (N). Fills the N column if absent.")
    nsz.add_argument("--n-cases",    type=int, default=None, metavar="INT",
                     help="Number of cases. Fills the N_cases column if absent.")
    nsz.add_argument("--n-controls", type=int, default=None, metavar="INT",
                     help="Number of controls. Fills the N_controls column if absent.")
    nsz.add_argument("--force-n",    action="store_true",
                     help="Overwrite existing N / N_cases / N_controls columns "
                          "with the values supplied via --n / --n-cases / --n-controls.")

    # ── QC thresholds ───────────────────────────────────────────────────────────
    qct = p.add_argument_group("QC thresholds (used when QC is enabled)")
    qct.add_argument("--eaf-min",  type=float, default=0.005, metavar="F",
                     help="EAF lower bound; symmetric upper = 1−F (default: 0.005).")
    qct.add_argument("--beta-max", type=float, default=5.0,   metavar="F",
                     help="Maximum |BETA| (default: 5.0).")
    qct.add_argument("--se-max",   type=float, default=5.0,   metavar="F",
                     help="Maximum SE (default: 5.0).")
    qct.add_argument("--info-min", type=float, default=0.4,   metavar="F",
                     help="Minimum INFO/imputation-quality score (default: 0.4).")
    qct.add_argument("--hwe-min",  type=float, default=1e-3,  metavar="F",
                     help="Minimum HWE p-value (default: 1e-3).")
    qct.add_argument("--mac-min",  type=int,   default=30,    metavar="N",
                     help="Minimum minor allele count (default: 30).")
    qct.add_argument("--daf-max",  type=float, default=0.12,  metavar="F",
                     help="Maximum |DAF|; set to 0 to skip DAF filter (default: 0.12).")

    return p.parse_args()

# ── Utrecht Science Park colour scheme ────────────────────────────────────────

UITHOF_COLOR = [
    "#FBB820", "#F59D10", "#E55738", "#DB003F", "#E35493", "#D5267B",
    "#CC0071", "#A8448A", "#9A3480", "#8D5B9A", "#705296", "#686AA9",
    "#6173AD", "#4C81BF", "#2F8BC9", "#1290D9", "#1396D8", "#15A6C1",
    "#5EB17F", "#86B833", "#C5D220", "#9FC228", "#78B113", "#49A01D",
    "#595A5C", "#A2A3A4", "#D7D8D7", "#ECECEC", "#FFFFFF", "#000000",
]

UITHOF_COLOR_LEGEND = {
    "yellow": "#FBB820", "gold": "#F59D10", "salmon": "#E55738",
    "darkpink": "#DB003F", "lightpink": "#E35493", "pink": "#D5267B",
    "hardpink": "#CC0071", "lightpurple": "#A8448A", "purple": "#9A3480",
    "lavendel": "#8D5B9A", "bluepurple": "#705296", "purpleblue": "#686AA9",
    "lightpurpleblue": "#6173AD", "seablue": "#4C81BF", "skyblue": "#2F8BC9",
    "azurblue": "#1290D9", "lightazurblue": "#1396D8", "greenblue": "#15A6C1",
    "seaweedgreen": "#5EB17F", "yellowgreen": "#86B833",
    "lightmossgreen": "#C5D220", "mossgreen": "#9FC228",
    "lightgreen": "#78B113", "green": "#49A01D",
    "grey": "#595A5C", "lightgrey": "#A2A3A4",
    "midgrey": "#D7D8D7", "verylightgrey": "#ECECEC",
    "white": "#FFFFFF", "black": "#000000",
}


# ── Helper functions ───────────────────────────────────────────────────────────

def ensure_dir(path: str) -> None:
    """Create directory (and parents) if it does not exist."""
    os.makedirs(path, exist_ok=True)


def save_tsv_gz(df: pd.DataFrame, path: str) -> None:
    """Write *df* as a tab-separated file, then gzip it in-place."""
    df.to_csv(path, sep="\t", index=False)
    with open(path, "rb") as f_in, gzip.open(path + ".gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(path)


def resolve_column(df: pd.DataFrame, aliases: list) -> str | None:
    """Return the first matching column name from *df*, case-insensitively."""
    df_cols_lower = {c.lower(): c for c in df.columns}
    for alias in aliases:
        if alias.lower() in df_cols_lower:
            return df_cols_lower[alias.lower()]
    return None


def reformat_output(df: pd.DataFrame) -> pd.DataFrame:
    """Rename GWASLab internal column names back to MetaGWASToolKit convention."""
    df = df.rename(columns={
        "POS":  "BP",
        "EA":   "EffectAllele",
        "NEA":  "OtherAllele",
        "INFO": "Info",
        "BETA": "Beta",
    })
    if "VariantID" not in df.columns:
        df["VariantID"] = df.get("SNPID", df.get("MarkerOriginal", ""))
    if "BetaMinor" not in df.columns and "Beta" in df.columns:
        df["BetaMinor"] = df["Beta"]
    if "MAF" not in df.columns and "EAF" in df.columns:
        df["MAF"] = df["EAF"]
    desired = [
        "VariantID", "MarkerOriginal", "rsID", "CHR", "BP", "Strand",
        "EffectAllele", "OtherAllele", "MinorAllele", "MajorAllele",
        "EAF", "MAF", "MAC", "HWE_P", "Info",
        "Beta", "BetaMinor", "SE", "P",
        "N", "N_cases", "N_controls", "Imputed", "DAF",
    ]
    df = df[[c for c in desired if c in df.columns]]
    for col in df.select_dtypes(include="category").columns:
        df[col] = df[col].cat.add_categories(["NA"])
    return df.fillna("NA")


def coerce_numeric_cols(data: pd.DataFrame, cols: list) -> None:
    """Coerce specified columns to numeric, setting non-convertible values to NaN."""
    for col in cols:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")


def plot_histogram(data: pd.DataFrame, col: str, label: str,
                   color: str, save_path: str) -> None:
    """Save a histogram of *col* to *save_path*."""
    plt.figure()
    sns.histplot(data=data, x=col, bins=25, kde=False,
                 stat="frequency", color=color)
    plt.title(f"Histogram of {label}")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def ref_vcf_path(ref_dir: str, population: str, build: str) -> str | None:
    """
    Return the expected 1KG reference VCF path for *population* and *build*,
    or None if no VCF exists for that build (e.g. hg18).

    Filename conventions:
      hg19 → <POP>.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz   (1KGP phase 3 v5)
      hg38 → <POP>.ALL.split_norm_af.1kg_30x.hg38.vcf.gz   (1KG 30x)
      hg18 → no 1KG VCF available; returns None
    """
    if build == "19":
        fname = f"{population}.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
    elif build == "38":
        fname = f"{population}.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
    else:
        return None
    return os.path.join(ref_dir, fname)


def dbsnp_vcf_path(ref_dir: str, build: str) -> str:
    """Return the expected dbSNP VCF path for *build*."""
    fname = "GCF_000001405.25.gz" if build == "19" else "GCF_000001405.40.gz"
    return os.path.join(ref_dir, fname)


# BGZF EOF block — every valid BGZF file must end with these 28 bytes.
_BGZF_EOF = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff"
    b"\x06\x00\x42\x43\x02\x00\x1b\x00"
    b"\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)

def check_bgzf(path: str) -> bool:
    """
    Return True if *path* ends with the canonical BGZF EOF block.
    A False result means the file is plain gzip (or truncated) and will be
    rejected by pysam / htslib with 'no BGZF EOF marker'.
    """
    try:
        with open(path, "rb") as fh:
            fh.seek(-28, 2)
            return fh.read(28) == _BGZF_EOF
    except (OSError, ValueError):
        return False


def assert_bgzf(path: str, label: str = "file") -> None:
    """
    Raise SystemExit with a helpful message if *path* is not BGZF-compressed.
    Call this before any step that passes the file to pysam / htslib.
    """
    if not check_bgzf(path):
        logging.error(
            "%s is not BGZF-compressed (pysam requires BGZF, not plain gzip):\n"
            "  %s\n"
            "Re-compress and index with:\n"
            "  gunzip -c '%s' | bgzip > '%s.bgzf.gz'\n"
            "  mv '%s.bgzf.gz' '%s'\n"
            "  tabix -p vcf '%s'",
            label, path, path, path, path, path, path,
        )
        sys.exit(1)


def build_qc_filter_expr(eaf: float, beta: float | None,
                          se: float | None, info: float | None,
                          daf: float | None) -> str:
    """Construct a gwaslab filter expression from QC thresholds."""
    filters = [f"(EAF >= {eaf} & EAF < {1 - eaf})"]
    if beta is not None:
        filters.append(f"(BETA >= {-beta} & BETA <= {beta})")
    if se is not None:
        filters.append(f"(SE <= {se})")
    if info is not None:
        filters.append(f"(INFO >= {info})")
    if daf is not None and daf > 0:
        filters.append(f"(DAF < {daf} & DAF > {-daf})")
    return " & ".join(filters)


def normalise_build(build: str) -> str:
    """Return canonical build string: '18', '19', or '38'."""
    return {
        "18": "18", "hg18": "18", "b36": "18", "grch36": "18", "36": "18",
        "19": "19", "grch37": "19", "hg19": "19", "b37": "19", "37": "19",
        "38": "38", "grch38": "38", "hg38": "38", "b38": "38",
    }.get(build.lower(), build.lower())


# ── Column alias tables ────────────────────────────────────────────────────────

SUMSTATS_ALIASES = {
    "snpid": ("SNPID", ["variantid", "snpid", "snp", "marker", "markername", "id", "variant"]),
    "chrom": ("CHR",   ["chr", "chrom", "chromosome", "chr(gcf1405.25)"]),
    "pos":   ("POS",   ["bp", "pos", "position", "base_pair_location", "bp_hg18", "bp_hg19",
                         "start(gcf1405.25)", "position_b38", "position_hg19", "position_hg38",
                         "bp_b37", "bp_b38", "pos_b37", "pos_b38"]),
    "ea":    ("EA",    ["effectallele", "ea", "a1", "allele1", "alt", "reference_allele",
                         "effect_allele", "riskallele", "codedallele"]),
    "nea":   ("NEA",   ["otherallele", "nea", "a2", "allele2", "ref", "non_effect_allele",
                         "other_allele", "noneffect_allele", "nonriskallele"]),
    "eaf":   ("EAF",   ["eaf", "effect_allele_frequency", "raf", "af", "allele_frequency",
                         "freq", "ref_allele_frequency", "effect_allele_freq", "caf",
                         "freq1", "freq(a1)", "freq.a1.1000g.eur", "a1_freq_1000g_eur"]),
    "beta":  ("BETA",  ["beta", "effect_size", "effectsize", "effect", "log_odds", "logor",
                         "beta_fixed", "b"]),
    "se":    ("SE",    ["se", "stderr", "standard_error", "sebeta", "log_odds_se",
                         "se_gc", "se_fixed"]),
    "p":     ("P",     ["p", "pval", "p_value", "pvalue", "p-value", "p-value_gc",
                         "p.value", "p_fixed"]),
    "n":     ("N",     ["n", "samplesize", "sample_size", "n_total", "ntotal", "n_samples",
                         "totalsamplesize", "n_eff", "neff"]),
    "rsid":  ("rsID",  ["rsid", "rs", "snp_id"]),
    "info":  ("INFO",  ["info", "impinfo", "imputation_quality", "r2", "rsq"]),
}

OPTIONAL_OTHER_ALIASES = {
    "CAVEAT":         ["caveat"],
    "HWE_P":          ["hwe_p", "hwep", "hwe"],
    "N_cases":        ["n_cases", "ncases", "cases", "n_case", "totalcases", "ncase",
                       "n_events", "n_event", "nevents", "nevent"],
    "N_controls":     ["n_controls", "ncontrols", "controls", "n_control", "ncontrol"],
    "MAF":            ["maf", "minor_allele_frequency", "minorallelefreq"],
    "MAC":            ["mac", "minor_allele_count"],
    "Strand":         ["strand"],
    "MarkerOriginal": ["markeroriginal", "marker_original"],
    "BetaMinor":      ["betaminor", "beta_minor"],
    "Imputed":        ["imputed"],
    "MajorAllele":    ["majorallele", "major_allele"],
    "MinorAllele":    ["minorallele", "minor_allele"],
    "DAF":            ["daf", "derived_allele_frequency"],
}

REQUIRED_COLS = {"snpid", "chrom", "pos", "ea", "nea", "beta", "se", "p"}


def file_tag(phenotype: str, population: str,
             input_build: str, output_build: str,
             added_n: bool = False) -> str:
    """
    Return the shared filename stem used by ALL output files and plots.

    Pattern:
      {phenotype}.{population}.input_b{input_build}.output_hg{output_build}.gwaslab
      …with optional trailing .added_n when N was supplied via --n / --n-cases / --n-controls

    Examples:
      LDL.EUR.input_b19.output_hg38.gwaslab
      CAD_SCHUNKERT.EUR.input_b18.output_hg38.gwaslab.added_n
      HDL.PAN.input_b38.output_hg38.gwaslab
    """
    stem = (f"{phenotype}.{population}"
            f".input_b{input_build}.output_hg{output_build}"
            f".gwaslab")
    if added_n:
        stem += ".added_n"
    return stem


# ── COJO output ────────────────────────────────────────────────────────────────

def write_cojo(gwas_obj, phenotype: str, population: str,
               input_build: str, output_build: str,
               output_loc: str, snpid_fmt: str, add_pos: bool,
               suffix: str = "", added_n: bool = False) -> None:
    """
    Write a GCTA-COJO format file.

    COJO spec columns:  SNP  A1  A2  freq  b  se  p  n
    With --cojo-pos:    SNP  CHR  BP  A1  A2  freq  b  se  p  n

    Parameters
    ----------
    snpid_fmt : 'chrpos' (CHR:BP) or 'rsid' (rsID column)
    add_pos   : prepend CHR and BP columns immediately after SNP
    suffix    : optional filename tag, e.g. 'qc'
    """
    df = gwas_obj.data   # read-only — no in-place modifications below

    # Build the SNP identifier column
    if snpid_fmt == "rsid":
        if "rsID" in df.columns:
            snp_col = df["rsID"].astype(str)
            logging.info("COJO SNP column: using rsID.")
        elif "SNPID" in df.columns:
            logging.warning(
                "COJO: --cojo-id rsid requested but no rsID column found; "
                "falling back to SNPID."
            )
            snp_col = df["SNPID"].astype(str)
        else:
            logging.warning(
                "COJO: --cojo-id rsid but no rsID or SNPID column; "
                "constructing CHR:BP instead."
            )
            snp_col = df["CHR"].astype(str) + ":" + df["POS"].astype(str)
    else:  # chrpos (default)
        snp_col = df["CHR"].astype(str) + ":" + df["POS"].astype(str)
        logging.info("COJO SNP column: using CHR:BP.")

    # Guard required COJO columns
    for req in ["EA", "NEA", "EAF", "BETA", "SE", "P", "N"]:
        if req not in df.columns:
            logging.error(
                "COJO output requires column '%s' which is absent. "
                "Skipping COJO file.", req
            )
            return

    cojo = pd.DataFrame({
        "SNP":  snp_col,
        "A1":   df["EA"],
        "A2":   df["NEA"],
        "freq": df["EAF"],
        "b":    df["BETA"],
        "se":   df["SE"],
        "p":    df["P"],
        "n":    df["N"],
    })

    if add_pos:
        # Insert CHR and BP right after the SNP column
        cojo.insert(1, "CHR", df["CHR"])
        cojo.insert(2, "BP",  df["POS"])

    tag       = f".{suffix}" if suffix else ""
    cojo_path = os.path.join(
        output_loc,
        f"{file_tag(phenotype, population, input_build, output_build, added_n)}{tag}.cojo.gz",
    )
    cojo.to_csv(cojo_path, sep="\t", index=False, compression="gzip")
    logging.info("[SAVE] COJO  → %s  (%s variants)", cojo_path, f"{len(cojo):,}")

# ── Stage checkpoint I/O ───────────────────────────────────────────────────────

def save_preprocess_checkpoint(gwas_data: pd.DataFrame, meta: dict,
                                stem: str, output_loc: str) -> None:
    """
    Save the standardised DataFrame + metadata JSON for the 'process' stage to load.

    Files written:
      {stem}.preprocess.parquet  — standardised DataFrame (BROTLI-compressed)
      {stem}.preprocess.json     — build metadata (reference, build_num, input_build)
    """
    parquet_path = os.path.join(output_loc, f"{stem}.preprocess.parquet")
    json_path    = os.path.join(output_loc, f"{stem}.preprocess.json")
    pq.write_table(pa.Table.from_pandas(gwas_data), parquet_path, compression="BROTLI")
    with open(json_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    logging.info("[SAVE] Preprocess parquet  → %s  (%s variants)",
                 parquet_path, f"{len(gwas_data):,}")
    logging.info("[SAVE] Preprocess metadata → %s", json_path)


def load_preprocess_checkpoint(stem: str, output_loc: str) -> tuple:
    """
    Load the standardised DataFrame + metadata JSON written by --stage preprocess.

    Returns (gwas_data, meta) where meta contains at minimum:
      reference, build_num, input_build
    """
    parquet_path = os.path.join(output_loc, f"{stem}.preprocess.parquet")
    json_path    = os.path.join(output_loc, f"{stem}.preprocess.json")
    for path, label in [(parquet_path, "Preprocess parquet"),
                        (json_path,    "Preprocess metadata")]:
        if not os.path.isfile(path):
            logging.error("%s checkpoint not found: %s\n"
                          "Run --stage preprocess first.", label, path)
            sys.exit(1)
    gwas_data = pq.read_table(parquet_path).to_pandas()
    with open(json_path) as fh:
        meta = json.load(fh)
    logging.info("[LOAD] Preprocess checkpoint: %s  (%s variants, build=%s)",
                 parquet_path, f"{len(gwas_data):,}", meta.get("reference", "?"))
    return gwas_data, meta


# ── Process sub-stage checkpoint I/O ──────────────────────────────────────────

# Suffix → (stage name for error messages, previous stage that writes it)
_PROCESS_CHECKPOINT_META = {
    "normalize":   ("process-normalize",    "preprocess"),
    "checkref":    ("process-check-ref",    "process-normalize"),
    "inferstrand": ("process-infer-strand", "process-check-ref"),
    "assignrsid":  ("process-assign-rsid",  "process-infer-strand"),
}


def save_process_checkpoint(gwas_obj, stem: str, output_loc: str, suffix: str) -> str:
    """
    Persist an intermediate Sumstats object as a pickle checkpoint.
    Returns the full path to the written pickle.

    File written: {stem}.{suffix}.pkl
    """
    pkl_path = os.path.join(output_loc, f"{stem}.{suffix}.pkl")
    logging.info("[SAVE] Process checkpoint (%s) → %s", suffix, pkl_path)
    gl.dump_pickle(gwas_obj, pkl_path, overwrite=True)
    return pkl_path


def load_process_checkpoint(stem: str, output_loc: str, suffix: str) -> object:
    """
    Load an intermediate Sumstats pickle written by a process sub-stage.
    Exits with an error message if the file is missing.
    """
    pkl_path = os.path.join(output_loc, f"{stem}.{suffix}.pkl")
    stage_name, prev_stage = _PROCESS_CHECKPOINT_META.get(
        suffix, (f"process-{suffix}", "prior stage")
    )
    if not os.path.isfile(pkl_path):
        logging.error(
            "[%s] Checkpoint not found: %s\n"
            "Run --stage %s first (with the same --gwas / --build / "
            "--liftover / --output / --dbsnp flags).",
            stage_name, pkl_path, prev_stage,
        )
        sys.exit(1)
    logging.info("[LOAD] Process checkpoint (%s): %s", suffix, pkl_path)
    gwas_obj = gl.load_pickle(pkl_path)
    if gwas_obj is None:
        logging.error("[%s] Failed to load checkpoint: %s", stage_name, pkl_path)
        sys.exit(1)
    return gwas_obj


# ── Per-chromosome checkpoint helpers ─────────────────────────────────────────

# gwaslab stores CHR as Int64 (nullable integer) after basic_check/_fix_chr:
#   1–22 = autosomes   23 = X   24 = Y   25 = nonPAR   26 = MT
# The SLURM array task ID maps directly to these integers, so no translation
# is needed — task 23 processes CHR == 23 (X chromosome).

def load_chrom_parquet(stem: str, output_loc: str, task_id: int,
                       suffix: str) -> pd.DataFrame:
    """
    Load a per-chromosome parquet written by split_by_chrom or a prior chr stage.

    If the file does not exist the chromosome is not present in this study;
    we log an info message and exit(0) so SLURM afterok dependencies are
    satisfied for the downstream array task.
    """
    path = os.path.join(output_loc, f"{stem}.chr{task_id}.{suffix}.parquet")
    if not os.path.isfile(path):
        logging.info(
            "[chr%d] Parquet not found — chromosome %d not present in this "
            "study. Exiting gracefully (exit 0).", task_id, task_id,
        )
        sys.exit(0)
    df = pq.read_table(path).to_pandas()
    logging.info("[chr%d] Loaded %s variants from %s", task_id, f"{len(df):,}", path)
    return df


def save_chrom_parquet(df: pd.DataFrame, stem: str, output_loc: str,
                       task_id: int, suffix: str) -> str:
    """
    Save a per-chromosome DataFrame as a BROTLI-compressed parquet.
    Returns the full path written.
    """
    path = os.path.join(output_loc, f"{stem}.chr{task_id}.{suffix}.parquet")
    pq.write_table(pa.Table.from_pandas(df), path, compression="BROTLI")
    logging.info("[chr%d] Saved %s variants → %s", task_id, f"{len(df):,}", path)
    return path


def make_sumstats_from_chrom_df(df: pd.DataFrame, reference: str) -> "gl.Sumstats":
    """
    Recreate a gwaslab Sumstats object from a per-chromosome DataFrame.

    Calls make_sumstats_object() then restores the STATUS column that was
    preserved through the parquet checkpoint.  STATUS is gwaslab's internal
    variant-state bitmask; if it is not carried over, gwaslab treats every
    variant as unprocessed and may re-run steps incorrectly.
    """
    status_backup = df["STATUS"].copy() if "STATUS" in df.columns else None
    gwas_obj = make_sumstats_object(df, reference)
    if status_backup is not None:
        gwas_obj.data["STATUS"] = status_backup.values
        logging.debug("STATUS column restored from parquet checkpoint (%d variants).",
                      len(status_backup))
    return gwas_obj


def split_by_chrom(gwas_obj, stem: str, output_loc: str) -> dict:
    """
    Split gwas_obj.data by CHR and save one parquet per chromosome.

    File names: {stem}.chr{N}.normalize.parquet  (N = CHR integer, 1–26)
    Manifest  : {stem}.chrsplit.json

    Returns the manifest dict (task_id_str → {chrom_int, n_variants}).
    """
    data = gwas_obj.data
    chr_values = data["CHR"].dropna().unique()
    logging.info("Splitting %s variants across %d chromosomes.",
                 f"{len(data):,}", len(chr_values))

    manifest: dict = {}
    for chrom_val in sorted(chr_values):
        task_id = int(chrom_val)
        subset  = data[data["CHR"] == chrom_val].reset_index(drop=True)
        if len(subset) == 0:
            continue
        save_chrom_parquet(subset, stem, output_loc, task_id, "normalize")
        manifest[str(task_id)] = {"chrom_int": task_id, "n_variants": len(subset)}

    json_path = os.path.join(output_loc, f"{stem}.chrsplit.json")
    with open(json_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    logging.info("[SAVE] Chr-split manifest → %s  (%d chromosomes)",
                 json_path, len(manifest))
    return manifest


def load_chrsplit_manifest(stem: str, output_loc: str) -> dict:
    """
    Load the chromosome manifest written by process-split.
    Exits with an error if the manifest is missing.
    """
    json_path = os.path.join(output_loc, f"{stem}.chrsplit.json")
    if not os.path.isfile(json_path):
        logging.error(
            "Chr-split manifest not found: %s\n"
            "Run --stage process-split first.", json_path,
        )
        sys.exit(1)
    with open(json_path) as fh:
        manifest = json.load(fh)
    logging.info("[LOAD] Chr-split manifest: %s  (%d chromosomes)",
                 json_path, len(manifest))
    return manifest


def run_merge(stem: str, output_loc: str, reference: str,
              build_num: str, input_build: str, args) -> None:
    """
    Merge stage: load all per-chromosome checkaf parquets, concatenate,
    recreate Sumstats (STATUS preserved), then run QC + plots + leads + COJO.

    This stage replaces the combination of process-check-af terminal outputs +
    the separate qc and cojo stages in the per-chromosome pipeline path.
    """
    manifest = load_chrsplit_manifest(stem, output_loc)

    # Choose which per-chr suffix to load: checkaf is always the terminal stage
    load_suffix = "checkaf"
    chr_dfs = []
    missing = []
    for task_id_str in sorted(manifest.keys(), key=int):
        task_id = int(task_id_str)
        path = os.path.join(output_loc, f"{stem}.chr{task_id}.{load_suffix}.parquet")
        if not os.path.isfile(path):
            missing.append(task_id)
            logging.warning("[merge] Missing chr %d parquet: %s", task_id, path)
            continue
        df = pq.read_table(path).to_pandas()
        chr_dfs.append(df)
        logging.info("[merge] Loaded chr %d: %s variants", task_id, f"{len(df):,}")

    if not chr_dfs:
        logging.error("[merge] No per-chromosome parquets found — cannot merge.")
        sys.exit(1)
    if missing:
        logging.warning("[merge] %d chromosome(s) missing from merge: %s",
                        len(missing), missing)

    combined = pd.concat(chr_dfs, ignore_index=True)
    logging.info("[merge] Combined: %s variants across %d chromosomes.",
                 f"{len(combined):,}", len(chr_dfs))

    added_n = any(v is not None for v in [args.n, args.n_cases, args.n_controls])
    plots_loc = os.path.join(output_loc, "PLOTS")
    ensure_dir(plots_loc)

    gwas_obj = make_sumstats_from_chrom_df(combined, build_num)
    del combined
    gc.collect()

    # Save raw (pre-QC) merged outputs
    save_raw_outputs(gwas_obj, args.gwas, args.population,
                     input_build, build_num, output_loc, added_n,
                     save_pickle=not args.no_pickle)

    if args.figures:
        stem_out = file_tag(args.gwas, args.population, input_build, build_num, added_n)
        plot_full_dataset(gwas_obj, args.gwas, build_num,
                          plots_loc, args.daf_max, stem_out)

    if args.cojo:
        write_cojo(gwas_obj, args.gwas, args.population,
                   input_build, build_num, output_loc,
                   snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
                   suffix="", added_n=added_n)

    # QC
    if args.qc:
        stem_out = file_tag(args.gwas, args.population, input_build, build_num, added_n)
        gwas_obj_qc = apply_qc(gwas_obj, args.eaf_min, args.beta_max,
                               args.se_max, args.info_min, args.daf_max)
        del gwas_obj
        gc.collect()

        save_qc_outputs(gwas_obj_qc, args.gwas, args.population,
                        input_build, build_num, output_loc, added_n,
                        save_pickle=not args.no_pickle)

        if args.figures:
            plot_qc_dataset(gwas_obj_qc, args.gwas, build_num,
                            plots_loc, args.daf_max, stem_out)

        if args.cojo:
            write_cojo(gwas_obj_qc, args.gwas, args.population,
                       input_build, build_num, output_loc,
                       snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
                       suffix="qc", added_n=added_n)

        if args.leads:
            extract_leads(gwas_obj_qc, args.gwas, args.population,
                          input_build, build_num, output_loc, True, added_n)
    else:
        if args.leads:
            extract_leads(gwas_obj, args.gwas, args.population,
                          input_build, build_num, output_loc, False, added_n)

    logging.info("[merge] Stage complete.")


# ── Process sub-stage runner functions ────────────────────────────────────────

def run_normalize(gwas_obj, reference: str, ref_loc: str,
                  n_cores: int, do_liftover: bool, population: str) -> tuple:
    """
    process-normalize: basic_check + remove_dup + liftover.
    Returns (gwas_obj, reference) — reference may be updated to '38' after liftover.
    """
    logging.info("\n===== Running basic_check =====")
    gwas_obj.basic_check(verbose=True)

    logging.info("\n===== Running remove_dup =====")
    gwas_obj.remove_dup(mode="md", keep_col="P", keep="first")
    logging.info("After duplicate removal: %d variants remain.", len(gwas_obj.data))

    logging.info("\n===== Running liftover =====")
    if do_liftover:
        ref_norm = normalise_build(reference)
        if ref_norm == "18":
            chain18 = os.path.join(ref_loc, "hg18ToHg38.over.chain.gz")
            if not os.path.isfile(chain18):
                logging.error(
                    "hg18→hg38 liftover requires 'hg18ToHg38.over.chain.gz' "
                    "in the reference directory:\n  %s", chain18,
                )
                sys.exit(1)
            logging.info("Lifting over from hg18 to hg38 using chain: %s …", chain18)
            gwas_obj.liftover(chain_path=chain18, to_build="38", remove=True)
            reference = "38"
            logging.info("Liftover complete. REFERENCE updated to '%s'.", reference)
        elif ref_norm == "19":
            logging.info("Lifting over from hg19 to hg38 …")
            gwas_obj.liftover(from_build="19", to_build="38", remove=True)
            reference = "38"
            logging.info("Liftover complete. REFERENCE updated to '%s'.", reference)
        elif ref_norm == "38":
            logging.info("Liftover skipped — data is already hg38.")
        else:
            logging.warning("Liftover skipped — unrecognised build '%s'.", reference)
    else:
        logging.info("Liftover skipped (--liftover not set).")

    return gwas_obj, reference


def run_check_ref(gwas_obj, reference: str, ref_loc: str) -> object:
    """
    process-check-ref: check_ref + flip_allele_stats + fix_id.
    """
    fasta_build = normalise_build(reference)
    fasta = os.path.join(ref_loc, f"hg{fasta_build}.fa.gz")

    logging.info("\n===== Running check_ref =====")
    if os.path.isfile(fasta):
        logging.info("Running check_ref with %s …", fasta)
        gwas_obj.check_ref(ref_seq=fasta)
        gwas_obj.flip_allele_stats()
    else:
        logging.warning("FASTA not found at '%s' — skipping check_ref.", fasta)

    logging.info("\n===== Running fix_id =====")
    gwas_obj.fix_id(fixid=True, forcefixid=True, overwrite=True)

    return gwas_obj


def run_infer_strand(gwas_obj, ref_loc: str, vcf: str, n_cores: int) -> object:
    """
    process-infer-strand: infer_strand2 + flip_allele_stats.
    High memory — bcftools sweep over the full 1KG VCF (~84 M variants).
    """
    logging.info("\n===== Running infer_strand2 =====")
    if vcf is None or not os.path.isfile(vcf):
        logging.warning("infer_strand2 skipped — reference VCF not available: %s", vcf)
    else:
        logging.info("Using reference file: %s …", vcf)
        gwas_obj.infer_strand2(vcf_path=vcf, threads=n_cores)
        gwas_obj.flip_allele_stats()

    return gwas_obj


def run_assign_rsid(gwas_obj, reference: str, ref_loc: str, n_cores: int) -> object:
    """
    process-assign-rsid: assign_rsid via dbSNP VCF sweep.
    Extreme memory — bcftools sweep over the full dbSNP VCF (~1 B variants for hg38).
    Only runs when --dbsnp is set; otherwise exits gracefully.
    """
    logging.info("\n===== Running assign_rsid =====")
    dbsnp = dbsnp_vcf_path(ref_loc, normalise_build(reference))
    if not os.path.isfile(dbsnp):
        logging.error("dbSNP VCF not found: %s — cannot run assign_rsid.", dbsnp)
        sys.exit(1)
    if not check_bgzf(dbsnp):
        logging.error(
            "dbSNP VCF is not BGZF-compressed (pysam requires BGZF, not plain gzip):\n"
            "  %s\nRe-compress with bgzip and index with tabix.", dbsnp,
        )
        sys.exit(1)
    logging.info("Using reference file: %s …", dbsnp)
    gwas_obj.harmonize(
        basic_check=False,
        ref_rsid_vcf=dbsnp,
        threads=n_cores,
        sweep_mode=True,
        verbose=True,
    )
    return gwas_obj


def run_check_af(gwas_obj, vcf: str, n_cores: int) -> object:
    """
    process-check-af: check_af2 (bcftools sweep over 1KG VCF).
    High memory — same sweep as infer_strand2 but against AF annotations.
    check_af2 is an annotation/flagging step; it does not flip alleles.
    """
    logging.info("\n===== Running check_af =====")
    if vcf is None or not os.path.isfile(vcf):
        logging.warning("check_af2 skipped — reference VCF not available: %s", vcf)
    else:
        gwas_obj.check_af2(vcf_path=vcf, ref_alt_freq="AF", threads=n_cores)

    return gwas_obj


# ── Pipeline steps ─────────────────────────────────────────────────────────────

# Detecting the column separator is a bit of a dark art, but we can make an 
# educated guess by looking at the first non-empty line of the file and 
# counting candidate delimiters.  This should work for most well-formed files, 
# and we can fall back to whitespace if nothing obvious is found.
def detect_separator(path: str) -> str:
    """
    Sniff the column separator from the first non-empty line of *path*.

    Reads the raw first line (handles .gz transparently) and counts
    candidate delimiters.  Returns the winner, falling back to whitespace
    (sep=r"\\s+") if nothing obvious is found.

    Candidates tested (in priority order):
        tab  →  "\t"
        pipe →  "|"
        semi →  ";"
        comma → ","
        space → r"\\s+"  (catches single-space and multi-space / mixed)
    """
    logging.info("\n===== Sniff out column separator =====")
    import gzip as _gzip

    opener = _gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if line:
                break
        else:
            return "\t"  # empty file — fall back

    candidates = [
        ("\t",      line.count("\t")),
        ("|",        line.count("|")),
        (";",        line.count(";")),
        (",",        line.count(",")),
    ]
    best_sep, best_n = max(candidates, key=lambda x: x[1])

    if best_n == 0:
        # No structured delimiter found — assume whitespace-separated
        logging.info("Separator sniff: no clear delimiter in header line; "
                     "using whitespace (sep=r'\\s+').")
        return r"\s+"

    sep_display = repr(best_sep)
    logging.info("Separator sniff: %s delimiter detected (%d occurrences in header).",
                 sep_display, best_n)
    return best_sep

# We can often infer the genome build from column names, e.g. "bp_hg19" or "position_b38".
# If such a hint is found, we compare it with the user-provided --build and log a warning 
# if they conflict. It also includes a sanity check on the POS column to flag suspiciously 
# low values that may indicate incorrect build or formatting issues.
def verify_build(gwas_data: pd.DataFrame, reference: str) -> str:
    """
    Infer build from column names and compare with *reference*.
    Returns the (possibly overridden) reference string.
    """
    logging.info("\n===== Verifying genome build =====")
    BUILD_HINTS_HG18 = [
        "bp_hg18", "position_hg18", "bp_b36", "position_b36",
        "bp_grch36", "position_grch36", "pos_b36",
    ]
    BUILD_HINTS_HG19 = [
        "chr(gcf1405.25)", "start(gcf1405.25)",
        "bp_hg19", "position_hg19",
        "bp_b37", "position_b37", "pos_b37", "bp_grch37", "position_grch37",
    ]
    BUILD_HINTS_HG38 = [
        "position_b38", "position_hg38", "bp_hg38",
        "bp_b38", "pos_b38", "bp_grch38", "position_grch38",
    ]

    ref_norm = {
        "18": "hg18", "hg18": "hg18", "b36": "hg18", "grch36": "hg18", "36": "hg18",
        "19": "hg19", "grch37": "hg19", "hg19": "hg19", "b37": "hg19", "37": "hg19",
        "38": "hg38", "grch38": "hg38", "hg38": "hg38", "b38": "hg38",
    }.get(reference.lower(), reference.lower())

    cols_lower = [c.lower() for c in gwas_data.columns]
    implied_build, hint_col = None, None

    for hint in BUILD_HINTS_HG18:
        if hint.lower() in cols_lower:
            implied_build, hint_col = "hg18", hint
            break
    if implied_build is None:
        for hint in BUILD_HINTS_HG19:
            if hint.lower() in cols_lower:
                implied_build, hint_col = "hg19", hint
                break
    if implied_build is None:
        for hint in BUILD_HINTS_HG38:
            if hint.lower() in cols_lower:
                implied_build, hint_col = "hg38", hint
                break

    if implied_build and implied_build != ref_norm:
        logging.warning("Column '%s' implies build '%s', but --build='%s' (%s). "
                        "Overriding to '%s'.", hint_col, implied_build, reference, ref_norm, implied_build)
        reference = implied_build
    elif implied_build:
        logging.info("Column '%s' confirms build '%s' — consistent with --build='%s'.",
                     hint_col, implied_build, reference)
    else:
        logging.info("No build-specific column names detected — assuming --build='%s'.", reference)

    # POS sanity check
    pos_col = resolve_column(gwas_data, [
        "bp", "pos", "position", "base_pair_location",
        "bp_hg18", "bp_hg19", "start(gcf1405.25)",
        "position_b38", "position_hg38", "position_hg19",
        "bp_b37", "bp_b38", "pos_b37", "pos_b38",
    ])
    if pos_col is not None:
        max_pos = pd.to_numeric(gwas_data[pos_col], errors="coerce").max()
        if max_pos < 1000:
            logging.warning("Maximum POS value is %.0f — suspiciously low, "
                            "check whether positions are in a non-standard format.", max_pos)
        else:
            logging.info("POS range check: max position = %s — looks reasonable.",
                         f"{max_pos:,.0f}")

    logging.info("Build verification complete. Using REFERENCE='%s'.", reference)
    return reference

# This is a potentially slow step, so we only do it if there are missing EAF values or 
# no EAF column at all.  We look up allele frequencies from the reference VCF using 
# pysam's TabixFile interface, matching on CHR, POS, EA, and NEA. 
# If the variant is found in the reference, we fill in the EAF; otherwise it remains NaN.
def check_and_fill_eaf(gwas_data: pd.DataFrame, ref_path: str) -> pd.DataFrame:
    """
    Check for an EAF column and fill missing values (or create the column)
    by looking up allele frequencies from a tabix-indexed VCF.
    """
    logging.info("\n===== Checking and filling EAF =====")
    import pysam
 
    EAF_ALIASES = [
        "eaf", "effect_allele_frequency", "raf", "af", "allele_frequency",
        "freq", "ref_allele_frequency", "effect_allele_freq", "caf",
        "freq1", "freq(a1)", "freq.a1.1000g.eur", "a1_freq_1000g_eur",
    ]
    eaf_col   = resolve_column(gwas_data, EAF_ALIASES)
    needs_eaf = (eaf_col is None) or (gwas_data[eaf_col].isna().sum() > 0)

    if not needs_eaf:
        logging.info("EAF column '%s' present and complete — no lookup needed.", eaf_col)
        return gwas_data

    chrom_col = resolve_column(gwas_data, ["chr", "chrom", "chromosome", "chr(gcf1405.25)"])
    pos_col   = resolve_column(gwas_data, ["bp", "pos", "position", "base_pair_location",
                                            "bp_hg18", "bp_hg19", "start(gcf1405.25)",
                                            "position_b38", "position_hg38", "position_hg19",
                                            "bp_b37", "bp_b38", "pos_b37", "pos_b38"])
    ea_col    = resolve_column(gwas_data, ["effectallele", "ea", "a1", "allele1", "alt",
                                            "reference_allele", "effect_allele", "riskallele", "codedallele"])
    nea_col   = resolve_column(gwas_data, ["otherallele", "nea", "a2", "allele2", "ref",
                                            "non_effect_allele", "other_allele", "noneffect_allele",
                                            "nonriskallele"])

    if any(c is None for c in [chrom_col, pos_col, ea_col, nea_col]):
        logging.warning("Cannot perform EAF lookup — missing required columns: "
                        "CHR=%s, POS=%s, EA=%s, NEA=%s. Skipping.", chrom_col, pos_col, ea_col, nea_col)
        return gwas_data

    logging.info("Opening reference VCF for EAF lookup: %s", ref_path)
    tbx = pysam.TabixFile(ref_path)

    def lookup_af(chrom, pos, ea, nea):
        """Look up AF for a single variant; returns NaN if not found."""
        try:
            for rec in tbx.fetch(str(chrom), int(pos) - 1, int(pos)):
                fields = rec.split("\t")
                ref, alt = fields[3], fields[4]
                if "," in alt:  # skip multi-allelic
                    continue
                info = dict(f.split("=") for f in fields[7].split(";") if "=" in f)
                af   = float(info.get("AF", "nan"))
                if alt == ea and ref == nea:
                    return af
                elif alt == nea and ref == ea:
                    return 1.0 - af
        except Exception:
            pass
        return np.nan

    if eaf_col is not None:
        n_missing = gwas_data[eaf_col].isna().sum()
        logging.info("Filling %s missing EAF values in '%s' from reference VCF.",
                     f"{n_missing:,}", eaf_col)
        mask = gwas_data[eaf_col].isna()
        gwas_data.loc[mask, eaf_col] = [
            lookup_af(r[chrom_col], r[pos_col], r[ea_col], r[nea_col])
            for _, r in gwas_data[mask].iterrows()
        ]
        n_still = gwas_data[eaf_col].isna().sum()
        logging.info("EAF filled for %s variants; %s still missing.",
                     f"{n_missing - n_still:,}", f"{n_still:,}")
    else:
        logging.info("Retrieving EAF for all %s variants from reference VCF.",
                     f"{len(gwas_data):,}")
        gwas_data["EAF"] = [
            lookup_af(r[chrom_col], r[pos_col], r[ea_col], r[nea_col])
            for _, r in gwas_data.iterrows()
        ]
        n_still = gwas_data["EAF"].isna().sum()
        logging.info("EAF retrieved: %s found, %s still missing.",
                     f"{gwas_data['EAF'].notna().sum():,}", f"{n_still:,}")

    tbx.close()
    return gwas_data

def apply_fixed_n(gwas_data: pd.DataFrame,
                  n: int | None, n_cases: int | None, n_controls: int | None,
                  force: bool) -> pd.DataFrame:
    """
    Fill N / N_cases / N_controls from command-line values.

    By default only fills columns that are absent or entirely NaN.
    With force=True, overwrites existing values.
    """
    pairs = [("N", n), ("N_cases", n_cases), ("N_controls", n_controls)]
    for col, val in pairs:
        if val is None:
            continue
        col_exists = col in gwas_data.columns and gwas_data[col].notna().any()
        if col_exists and not force:
            logging.info("Fixed N: column '%s' already present — skipping "
                         "(use --force-n to override).", col)
        else:
            if col_exists:
                logging.warning("Fixed N: overwriting existing '%s' column with %d "
                                "(--force-n set).", col, val)
            else:
                logging.info("Fixed N: setting '%s' = %d for all variants.", col, val)
            gwas_data[col] = val
    return gwas_data


# This function checks for the presence of key columns and attempts to correct or
# derive them if possible.
def correct_columns(gwas_data: pd.DataFrame) -> pd.DataFrame:
    """Correct / derive CAVEAT, P, SE, and N columns."""
    logging.info("\n===== Column Correction =====")

    # CAVEAT
    if "CAVEAT" not in gwas_data.columns:
        gwas_data["CAVEAT"] = "None"
    gwas_data["CAVEAT"] = gwas_data["CAVEAT"].fillna("None")
    logging.info("CAVEAT column ensured.")

    # log(P) → P
    p_col    = resolve_column(gwas_data, ["p", "pval", "p_value", "pvalue",
                                           "p-value", "p-value_gc", "p.value", "p_fixed"])
    logp_col = resolve_column(gwas_data, ["log_p", "logp", "log10p", "log10_p", "mlog10p",
                                           "log_pvalue", "log_p_value", "-log10p",
                                           "neg_log10_p", "neg_log_p", "p_log", "log(p)"])
    if p_col is None and logp_col is not None:
        col_values = gwas_data[logp_col].astype(float)
        if (col_values < 0).sum() > (col_values > 0).sum():
            logging.info("'%s' appears to be log10(P) (mostly negative) — converting via 10^x.", logp_col)
            gwas_data["P"] = 10 ** col_values
        else:
            logging.info("'%s' appears to be -log10(P) (mostly positive) — converting via 10^(-x).", logp_col)
            gwas_data["P"] = 10 ** (-col_values)
        gwas_data["P"] = gwas_data["P"].clip(lower=5e-324, upper=1.0)
        n_invalid = ((gwas_data["P"] <= 0) | (gwas_data["P"] > 1)).sum()
        if n_invalid:
            logging.warning("%d P-values outside (0, 1] after conversion.", n_invalid)
        else:
            logging.info("P converted to column 'P' — all values valid.")
    elif p_col is not None:
        logging.info("P-value column found: '%s' — no conversion needed.", p_col)
    else:
        logging.warning("No P or log(P) column found.")

    # SE back-calculation
    se_col   = resolve_column(gwas_data, ["se", "stderr", "standard_error", "sebeta",
                                           "log_odds_se", "se_gc", "se_fixed"])
    beta_col = resolve_column(gwas_data, ["beta", "effect_size", "effectsize", "effect",
                                           "log_odds", "logor", "beta_fixed", "b"])
    p_col    = resolve_column(gwas_data, ["p", "pval", "p_value", "pvalue",
                                           "p-value", "p-value_gc", "p.value", "p_fixed"])
    if se_col is not None:
        logging.info("SE column found: '%s' — no back-calculation needed.", se_col)
    elif beta_col is not None and p_col is not None:
        logging.info("SE not found — back-calculating from '%s' and '%s'.", beta_col, p_col)
        z = np.abs(norm.ppf(gwas_data[p_col].clip(lower=1e-300) / 2))
        z = np.where(z == 0, np.nan, z)
        gwas_data["SE"] = np.abs(gwas_data[beta_col]) / z
        logging.info("SE back-calculated and stored in column 'SE'.")
    else:
        logging.warning("SE not found and cannot be back-calculated — Beta and/or P column missing.")

    # N derivation
    n_col        = resolve_column(gwas_data, ["n", "samplesize", "sample_size", "n_total",
                                               "ntotal", "n_samples", "totalsamplesize", "n_eff", "neff"])
    ncase_col    = resolve_column(gwas_data, ["n_cases", "ncases", "cases", "n_case",
                                               "totalcases", "ncase", "n_events", "n_event",
                                               "nevents", "nevent"])
    ncontrol_col = resolve_column(gwas_data, ["n_controls", "ncontrols", "controls",
                                               "n_control", "ncontrol"])

    if n_col is not None and ncase_col is not None and ncontrol_col is not None:
        logging.info("N, N_cases, and N_controls all present — no derivation needed.")
    elif n_col is None and ncase_col is not None and ncontrol_col is not None:
        gwas_data["N"] = gwas_data[ncase_col] + gwas_data[ncontrol_col]
        logging.info("N calculated (median=%.0f) from '%s' + '%s'.",
                     gwas_data["N"].median(), ncase_col, ncontrol_col)
    elif n_col is not None and ncase_col is not None and ncontrol_col is None:
        gwas_data["N_controls"] = gwas_data[n_col] - gwas_data[ncase_col]
        logging.info("N_controls calculated (median=%.0f) from '%s' - '%s'.",
                     gwas_data["N_controls"].median(), n_col, ncase_col)
    elif n_col is not None and ncontrol_col is not None and ncase_col is None:
        gwas_data["N_cases"] = gwas_data[n_col] - gwas_data[ncontrol_col]
        logging.info("N_cases calculated (median=%.0f) from '%s' - '%s'.",
                     gwas_data["N_cases"].median(), n_col, ncontrol_col)
    else:
        logging.warning("Cannot derive N — insufficient sample size columns available "
                        "(N=%s, N_cases=%s, N_controls=%s).", n_col, ncase_col, ncontrol_col)

    logging.info("Column correction complete. Shape: %s", gwas_data.shape)
    return gwas_data

# This function performs the core column standardisation by matching against the alias tables.
def standardise_columns(gwas_data: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to GWASLab standard names."""
    logging.info("\n===== Column Standardisation =====")

    rename_map, missing_required = {}, []

    for canonical, (standard_name, aliases) in SUMSTATS_ALIASES.items():
        matched = resolve_column(gwas_data, aliases)
        if matched:
            if matched != standard_name:
                rename_map[matched] = standard_name
                logging.info("Renaming '%s' -> '%s'.", matched, standard_name)
            else:
                logging.info("'%s' already correctly named — no rename needed.", standard_name)
        elif canonical in REQUIRED_COLS:
            missing_required.append(canonical)
            logging.warning("Required column '%s' (-> '%s') not found — tried: %s.",
                            canonical, standard_name, aliases)
        else:
            logging.debug("Optional column '%s' not found — skipping.", canonical)

    for standard_name, aliases in OPTIONAL_OTHER_ALIASES.items():
        matched = resolve_column(gwas_data, aliases)
        if matched and matched != standard_name:
            rename_map[matched] = standard_name
            logging.info("Renaming '%s' -> '%s'.", matched, standard_name)

    if rename_map:
        gwas_data.rename(columns=rename_map, inplace=True)
        logging.info("Renamed %d column(s).", len(rename_map))
    else:
        logging.info("No columns needed renaming.")

    if missing_required:
        logging.warning("Missing required columns: %s. Sumstats creation may fail.", missing_required)

    logging.info("Standardisation complete. Columns now: %s", list(gwas_data.columns))
    return gwas_data


def plot_raw_histograms(gwas_data: pd.DataFrame, stem: str,
                        plots_loc: str) -> None:
    """Generate pre-processing histograms for key columns."""
    FIGURE_COLUMNS = {
        "BETA":       ("Effect Size (Beta)",      "#E55738"),
        "EAF":        ("Effect Allele Frequency", "#1290D9"),
        "SE":         ("Standard Error (SE)",     "#49A01D"),
        "INFO":       ("Imputation INFO Score",   "#705296"),
        "N":          ("Sample Size (N)",         "#F59D10"),
        "N_cases":    ("Number of Cases",         "#E35493"),
        "N_controls": ("Number of Controls",      "#595A5C"),
    }
    logging.info("\n===== Generating pre-processing histograms =====")
    for col, (label, color) in FIGURE_COLUMNS.items():
        if col in gwas_data.columns:
            logging.info("Plotting '%s' (%s).", col, label)
            plot_histogram(
                gwas_data, col, label, color,
                os.path.join(plots_loc, f"{stem}.histogram.{col}.raw.png"),
            )
        else:
            logging.info("Column '%s' not present — skipping histogram.", col)


def make_sumstats_object(gwas_data: pd.DataFrame, reference: str) -> "gl.Sumstats":
    """Create a gwaslab Sumstats object from *gwas_data*."""
    logging.info("\n===== Creating Sumstats Object =====")

    other_cols = [c for c in OPTIONAL_OTHER_ALIASES if c in gwas_data.columns]
    sumstats_kwargs = {
        canonical: standard_name
        for canonical, (standard_name, _) in SUMSTATS_ALIASES.items()
        if standard_name in gwas_data.columns
    }

    # EAF placeholder required by harmonize / infer_strand
    if "EAF" not in gwas_data.columns:
        logging.info("Adding EAF placeholder column (all NaN) for gwaslab compatibility.")
        gwas_data["EAF"] = float("nan")

    gwas_obj = gl.Sumstats(
        gwas_data,
        other=other_cols,
        build=reference,
        verbose=True,
        **sumstats_kwargs,
    )
    logging.info("Sumstats object created with %d variants.", len(gwas_obj.data))
    return gwas_obj

# This function runs the core processing steps in sequence, with logging and error handling.
def run_processing(gwas_obj, reference: str, ref_loc: str, vcf: str,
                   n_cores: int, do_liftover: bool, do_dbsnp: bool,
                   population: str = "EUR") -> tuple:
    """
    Run the core GWASLab processing steps.
    Returns (gwas_obj, reference) — reference may be updated after liftover.
    vcf is updated internally after liftover so downstream steps use the
    correct build's reference VCF.
    """
    # basic_check
    logging.info("\n===== Running basic_check =====")
    gwas_obj.basic_check(verbose=True)

    # remove_dup
    logging.info("\n===== Running remove_dup =====")
    gwas_obj.remove_dup(mode="md", keep_col="P", keep="first")
    logging.info("After duplicate removal: %d variants remain.", len(gwas_obj.data))

    # liftover
    logging.info("\n===== Running liftover =====")
    if do_liftover:
        ref_norm = normalise_build(reference)
        if ref_norm == "18":
            # hg18→hg38: gwaslab only ships hg19↔hg38 chains, so we require a
            # hg18ToHg38.over.chain.gz in the ref directory. A two-step approach
            # (18→19→38) is blocked by STATUS code tracking after the first pass.
            chain18 = os.path.join(ref_loc, "hg18ToHg38.over.chain.gz")
            if not os.path.isfile(chain18):
                logging.error(
                    "hg18→hg38 liftover requires 'hg18ToHg38.over.chain.gz' "
                    "in the reference directory:\n  %s\n"
                    "Download it from UCSC:\n"
                    "  wget -P '%s' https://hgdownload.soe.ucsc.edu/goldenPath/"
                    "hg18/liftOver/hg18ToHg38.over.chain.gz",
                    chain18, ref_loc,
                )
                sys.exit(1)
            logging.info("Lifting over from hg18 to hg38 using chain: %s …", chain18)
            gwas_obj.liftover(chain_path=chain18, to_build="38", remove=True)
            reference = "38"
            logging.info("Liftover complete. REFERENCE updated to '%s'.", reference)
        elif ref_norm == "19":
            logging.info("Lifting over from hg19 to hg38 …")
            gwas_obj.liftover(from_build="19", to_build="38", remove=True)
            reference = "38"
            logging.info("Liftover complete. REFERENCE updated to '%s'.", reference)
        elif ref_norm == "38":
            logging.info("Liftover skipped — data is already hg38.")
        else:
            logging.warning("Liftover skipped — unrecognised build '%s'.", reference)

        # Update VCF to the post-liftover build so infer_strand2 / check_af2
        # use the correct reference (hg38 after any liftover).
        # This also resolves the VCF for build 18 (which starts as None).
        new_vcf = ref_vcf_path(ref_loc, population, normalise_build(reference))
        if new_vcf is not None and os.path.isfile(new_vcf):
            if new_vcf != vcf:
                logging.info("Reference VCF updated to post-liftover build: %s", new_vcf)
            vcf = new_vcf
        else:
            logging.warning(
                "Post-liftover reference VCF not found: %s — "
                "continuing with pre-liftover VCF: %s", new_vcf, vcf,
            )
    else:
        logging.info("Liftover skipped (--liftover not set).")

    # check_ref + flip
    # NOTE: flip_allele_stats() is still required — check_ref (like
    # infer_strand/infer_strand2) only updates STATUS codes; the 
    # actual re-orientation of BETA/EAF/OR etc. is done by flip_allele_stats().
    logging.info("\n===== Running check_ref =====")
    fasta_build = normalise_build(reference)
    fasta = os.path.join(ref_loc, f"hg{fasta_build}.fa.gz")
    if os.path.isfile(fasta):
        logging.info("Running check_ref with %s …", fasta)
        gwas_obj.check_ref(ref_seq=fasta)
        # gc.collect() # for debug of killed 9 signal
        gwas_obj.flip_allele_stats()
    else:
        logging.warning("FASTA not found at '%s' — skipping check_ref.", fasta)

    # fix_id
    logging.info("\n===== Running fix_id =====")
    gwas_obj.fix_id(fixid=True, forcefixid=True, overwrite=True)

    # infer_strand2 + flip
    # infer_strand2 replaces infer_strand: instead of per-variant tabix queries
    # it uses bcftools to sweep the entire reference VCF once, building a lookup
    # table in memory — much faster and lower peak I/O for large datasets.
    # convert_to_bcf=True lets bcftools work on a BCF internally, which is
    # faster than repeated VCF.gz decompression.
    # NOTE: flip_allele_stats() is still required — infer_strand2 (like
    # infer_strand) only updates STATUS codes; the actual re-orientation of
    # BETA/EAF/OR etc. is done by flip_allele_stats().
    logging.info("\n===== Running infer_strand2 =====")
    if vcf is None:
        logging.warning("infer_strand2 skipped — no reference VCF available.")
    else:
        logging.info("Using reference file: %s …", vcf)
        gwas_obj.infer_strand2(vcf_path=vcf, threads=n_cores)
        gwas_obj.flip_allele_stats()

    # assign_rsid
    # Use harmonize(basic_check=False, sweep_mode=True) for the rsID-only step:
    # this routes through _assign_rsid() (bcftools one-pass sweep) instead of
    # _parallelize_assign_rsid() (per-chunk tabix). Substantially faster and
    # lower-memory for large datasets. All other harmonize steps are disabled
    # (basic_check=False, ref_seq/ref_infer=None) so only assign_rsid runs.
    logging.info("\n===== Running assign_rsid =====")
    if do_dbsnp:
        dbsnp = dbsnp_vcf_path(ref_loc, normalise_build(reference))
        if not os.path.isfile(dbsnp):
            logging.error("dbSNP VCF not found: %s — skipping rsID assignment.", dbsnp)
        elif not check_bgzf(dbsnp):
            logging.error(
                "dbSNP VCF is not BGZF-compressed (pysam requires BGZF, not plain gzip):\n"
                "  %s\n"
                "Re-compress and index with:\n"
                "  gunzip -c '%s' | bgzip > '%s.tmp.gz' && mv '%s.tmp.gz' '%s'\n"
                "  tabix -p vcf '%s'",
                dbsnp, dbsnp, dbsnp, dbsnp, dbsnp, dbsnp,
            )
            sys.exit(1)
        else:
            logging.info("Using reference file: %s …", dbsnp)
            gwas_obj.harmonize(
                basic_check=False,
                ref_rsid_vcf=dbsnp,
                threads=n_cores,
                sweep_mode=True,
                verbose=True,
            )
    else:
        logging.info("dbSNP rsID assignment skipped (--dbsnp not set).")

    # check_af2 (sweep variant of check_af)
    # Uses _check_af_with_annotation: one bcftools pass over the reference VCF
    # rather than per-chunk tabix queries (_parallelize_check_af).
    # check_af2 takes vcf_path= (not ref_infer=) and has no sweep_mode flag —
    # it is always the sweep path.
    logging.info("\n===== Running check_af =====")
    if vcf is None:
        logging.warning("check_af2 skipped — no reference VCF available.")
    else:
        gwas_obj.check_af2(vcf_path=vcf, ref_alt_freq="AF", threads=n_cores)

    return gwas_obj, reference

# This function saves the outputs at each stage in multiple formats, with logging.
def save_raw_outputs(gwas_obj, phenotype: str, population: str,
                     input_build: str, output_build: str,
                     output_loc: str, added_n: bool = False,
                     save_pickle: bool = True) -> str:
    """Save raw (pre-QC) outputs: pickle (optional), parquet, TSV.GZ. Returns parquet path."""
    stem         = file_tag(phenotype, population, input_build, output_build, added_n)
    pkl_path     = os.path.join(output_loc, f"{stem}.pkl")
    parquet_path = os.path.join(output_loc, f"{stem}.parquet")
    tsv_path     = os.path.join(output_loc, f"{stem}.tsv")

    logging.info("\n===== Saving raw (pre-QC) outputs =====")
    if save_pickle:
        logging.info("[SAVE] Pickle  → %s", pkl_path)
        gl.dump_pickle(gwas_obj, pkl_path, overwrite=True)
    else:
        logging.info("[SAVE] Pickle  skipped (--no-pickle).")
    gwas_obj.log.show()
    gwas_obj.log.save(pkl_path.replace(".pkl", ".log"))

    logging.info("[SAVE] Parquet → %s", parquet_path)
    pq.write_table(pa.Table.from_pandas(gwas_obj.data), parquet_path, compression="BROTLI")

    # Write TSV directly from gwas_obj.data — avoids reading the parquet back
    # into memory as a second copy just to reformat and write it out.
    logging.info("[SAVE] TSV.GZ  → %s.gz", tsv_path)
    save_tsv_gz(reformat_output(gwas_obj.data), tsv_path)

    return parquet_path

# This function generates the standard set of plots for the dataset before QC filtering.
def plot_full_dataset(gwas_obj, phenotype: str, reference: str,
                      plots_loc: str, daf_max: float, stem: str) -> None:
    """Generate Manhattan, QQ, and DAF plots for the full (pre-QC) dataset."""
    logging.info("\n===== Generating plots for the full (pre-QC) dataset =====")
    daf_thr = daf_max if daf_max > 0 else 0.12
    gwas_obj.plot_daf(
        threshold=daf_thr,
        save=os.path.join(plots_loc, f"{stem}.EAF.png"),
        save_kwargs={"dpi": 300},
    )
    for mode, suffix in [("m", "manhattan"), ("qq", "qq")]:
        # anno="GENENAME" triggers gwaslab's to_annotate code path, which is
        # only initialised inside the Manhattan/regional block.  Passing anno
        # on a pure QQ call causes UnboundLocalError in gwaslab ≤0.4.x.
        is_manhattan = mode == "m"
        gwas_obj.plot_mqq(
            skip=2, cut=10, mode=mode,
            sig_line=True, sig_level=5e-8,
            anno="GENENAME" if is_manhattan else None,
            anno_style="right" if is_manhattan else None,
            windowsizekb=500 if is_manhattan else None,
            arm_offset=2 if is_manhattan else None,
            repel_force=0.02 if is_manhattan else None,
            use_rank=True, build=reference,
            stratified=True, drop_chr_start=True,
            title=phenotype,
            save=os.path.join(plots_loc, f"{stem}.{suffix}.500kb.300dpi.png"),
            save_kwargs={"dpi": 300},
            verbose=True,
        )


def apply_qc(gwas_obj, eaf_min: float, beta_max: float, se_max: float,
             info_min: float, daf_max: float) -> "gl.Sumstats":
    """Apply numeric QC filters and return a filtered Sumstats object."""
    logging.info("\n===== Applying QC filters =====")
    cols = set(gwas_obj.data.columns)
    coerce_numeric_cols(gwas_obj.data,
                        ["EAF", "DAF", "BETA", "SE", "INFO", "MAC", "HWE_P",
                         "N", "N_cases", "N_controls", "POS"])
    if "DAF" in cols and daf_max > 0:
        gwas_obj.data["DAF"] = gwas_obj.data["DAF"].fillna(0.0)

    # Only include a filter term when the column is actually present in the
    # data — pandas query() raises UndefinedVariableError for absent columns.
    optional_filters = [
        ("BETA", beta_max,  f"--beta-max {beta_max}"),
        ("SE",   se_max,    f"--se-max {se_max}"),
        ("INFO", info_min,  f"--info-min {info_min}"),
        ("DAF",  daf_max,   f"--daf-max {daf_max}"),
    ]
    for col, val, flag in optional_filters:
        if col not in cols:
            logging.warning(
                "QC: column '%s' not present in data — skipping %s filter.", col, flag,
            )

    expr = build_qc_filter_expr(
        eaf_min,
        beta_max  if "BETA" in cols else None,
        se_max    if "SE"   in cols else None,
        info_min  if "INFO" in cols else None,
        daf_max   if "DAF"  in cols else 0,
    )
    logging.info("QC filter expression: %s", expr)
    gwas_obj_qc = gwas_obj.filter_value(expr=expr)
    logging.info("Variants after QC: %d", len(gwas_obj_qc.data))
    return gwas_obj_qc

# This function saves the QC-filtered outputs in multiple formats, with logging.
def save_qc_outputs(gwas_obj_qc, phenotype: str, population: str,
                    input_build: str, output_build: str,
                    output_loc: str, added_n: bool = False,
                    save_pickle: bool = True) -> str:
    """Save QC-filtered outputs: pickle (optional), parquet, TSV.GZ. Returns parquet path."""
    logging.info("\n===== Saving QC-filtered outputs =====")
    stem         = file_tag(phenotype, population, input_build, output_build, added_n)
    pkl_path     = os.path.join(output_loc, f"{stem}.qc.pkl")
    parquet_path = os.path.join(output_loc, f"{stem}.qc.parquet")
    tsv_path     = os.path.join(output_loc, f"{stem}.qc.tsv")

    if save_pickle:
        logging.info("[SAVE] QC Pickle  → %s", pkl_path)
        gl.dump_pickle(gwas_obj_qc, pkl_path, overwrite=True)
    else:
        logging.info("[SAVE] QC Pickle  skipped (--no-pickle).")
    gwas_obj_qc.log.save(pkl_path.replace(".pkl", ".log"))

    logging.info("[SAVE] QC Parquet → %s", parquet_path)
    pq.write_table(pa.Table.from_pandas(gwas_obj_qc.data), parquet_path, compression="BROTLI")

    # Write TSV directly from gwas_obj_qc.data — avoids reading the parquet back
    # into memory as a second copy just to reformat and write it out.
    logging.info("[SAVE] QC TSV.GZ  → %s.gz", tsv_path)
    save_tsv_gz(reformat_output(gwas_obj_qc.data), tsv_path)

    return parquet_path

# This function generates the standard set of plots for the QC-filtered dataset.
def plot_qc_dataset(gwas_obj_qc, phenotype: str, reference: str,
                    plots_loc: str, daf_max: float, stem: str) -> None:
    """Generate Manhattan, QQ, and DAF plots for the QC-filtered dataset."""
    logging.info("\n===== Generating plots for the QC-filtered dataset =====")
    daf_thr = daf_max if daf_max > 0 else 0.12
    gwas_obj_qc.plot_daf(
        threshold=daf_thr,
        save=os.path.join(plots_loc, f"{stem}.EAF.qc.png"),
        save_kwargs={"dpi": 300},
    )
    for mode, suffix in [("m", "manhattan"), ("qq", "qq")]:
        is_manhattan = mode == "m"
        gwas_obj_qc.plot_mqq(
            skip=2, cut=10, mode=mode,
            sig_line=True, sig_level=5e-8,
            anno="GENENAME" if is_manhattan else None,
            anno_style="right" if is_manhattan else None,
            windowsizekb=500 if is_manhattan else None,
            arm_offset=2 if is_manhattan else None,
            repel_force=0.02 if is_manhattan else None,
            use_rank=True, build=reference,
            stratified=True, drop_chr_start=True,
            title=f"{phenotype} (QC)",
            save=os.path.join(plots_loc, f"{stem}.{suffix}.500kb.300dpi.qc.png"),
            save_kwargs={"dpi": 300},
            verbose=True,
        )

# This function extracts the genome-wide significant lead variants from the 
# (QC-filtered) dataset and saves them as TSV.
def extract_leads(source_obj, phenotype: str, population: str,
                  input_build: str, output_build: str,
                  output_loc: str, is_qc: bool,
                  added_n: bool = False) -> None:
    """Extract and save genome-wide significant lead variants."""
    logging.info("\n===== Extracting lead variants =====")
    suffix    = "qc.leads" if is_qc else "leads"
    stem      = file_tag(phenotype, population, input_build, output_build, added_n)
    leads_tsv = os.path.join(output_loc, f"{stem}.{suffix}.tsv")

    leads = source_obj.get_lead(anno=True, build=output_build, sig_level=5e-8, verbose=True)

    if hasattr(leads, "to_csv"):
        leads.to_csv(leads_tsv, index=False)
    elif hasattr(leads, "data"):
        leads.data.to_csv(leads_tsv, index=False)

    logging.info("[SAVE] Lead variants → %s", leads_tsv)


# ── Main ───────────────────────────────────────────────────────────────────────
# The main function orchestrates the entire workflow.  It supports both the
# classic end-to-end run (--stage all, the default) and a staged execution
# model (--stage preprocess / process / qc / cojo) where each stage saves a
# checkpoint so successive stages can be submitted as separate SLURM jobs with
# individual memory and time limits.
#
# Stage handoff files (checkpoint chain):
#   preprocess          → process-normalize   : {stem}.preprocess.parquet + .json
#   process-normalize   → process-check-ref   : {stem}.normalize.pkl
#   process-check-ref   → process-infer-strand: {stem}.checkref.pkl
#   process-infer-strand→ process-assign-rsid  : {stem}.inferstrand.pkl
#   process-assign-rsid → process-check-af    : {stem}.assignrsid.pkl
#   process-check-af    → qc / cojo           : {stem}.pkl + .parquet + .tsv.gz
#   qc                  → cojo                : {stem}.qc.pkl + .qc.parquet + .qc.tsv.gz
#
# Pass identical --gwas / --build / --liftover / --output / --dbsnp flags to
# every stage so that file stems and checkpoint paths are consistent.
def main() -> None:
    args = parse_args()

    # ── Resolve paths ──────────────────────────────────────────────────────────
    input_path = (args.input if os.path.isabs(args.input)
                  else os.path.join(args.directory, args.input))

    if args.output:
        if "{" in args.output or "}" in args.output:
            print(
                f"[ERROR] --output contains unresolved placeholder: {args.output!r}\n"
                "        Do not pass literal brace expressions; omit --output to use\n"
                f"        the default: <--dir>/{args.gwas}/GWASCatalog",
                file=sys.stderr,
            )
            sys.exit(1)
        output_loc = args.output
    else:
        output_loc = os.path.join(args.directory, args.gwas, "GWASCatalog")

    plots_loc = os.path.join(output_loc, "PLOTS")
    ensure_dir(output_loc)
    ensure_dir(plots_loc)

    # ── Logging ────────────────────────────────────────────────────────────────
    log_path = os.path.join(output_loc, f"{args.gwas}.gwaslab_process.log")
    setup_logging(log_path)

    logging.info("=" * 60)
    logging.info("%s  v%s  (%s)", VERSION_NAME, VERSION, VERSION_DATE)
    logging.info("=" * 60)
    logging.info("Study        : %s", args.gwas)
    logging.info("Input file   : %s", input_path)
    logging.info("Reference dir: %s", args.ref)
    logging.info("Output dir   : %s", output_loc)
    logging.info("Population   : %s", args.population)
    logging.info("Build        : %s", args.build)
    logging.info("Toggles      : liftover=%s  dbsnp=%s  qc=%s  "
                 "only_qc=%s  fill_eaf=%s  figures=%s  leads=%s  "
                 "no_pickle=%s  stage=%s  chrom=%s  threads=%d",
                 args.liftover, args.dbsnp, args.qc,
                 args.only_qc, args.fill_eaf, args.figures, args.leads,
                 args.no_pickle, args.stage,
                 args.chrom if args.chrom else "(whole-genome)",
                 args.threads)
    if args.cojo:
        logging.info("COJO options : id=%s  pos=%s", args.cojo_id, args.cojo_pos)
    if any(v is not None for v in [args.n, args.n_cases, args.n_controls]):
        logging.info("Fixed N      : N=%s  N_cases=%s  N_controls=%s  force=%s",
                     args.n, args.n_cases, args.n_controls, args.force_n)
    logging.info("Threads      : %d", args.threads)

    if args.cojo and args.cojo_id == "rsid" and not args.dbsnp:
        logging.warning(
            "--cojo-id rsid selected but --dbsnp is not set. "
            "COJO SNP column will use rsID only if the input already contains one."
        )

    REFERENCE   = args.build
    input_build = normalise_build(args.build)   # fixed; never changes after liftover
    added_n     = any(v is not None for v in [args.n, args.n_cases, args.n_controls])

    # ── GWASLab data directory ─────────────────────────────────────────────────
    gl.options.set_option("data_directory", args.ref)

    # ── Resolve reference VCF ──────────────────────────────────────────────────
    build_num = normalise_build(REFERENCE)
    vcf = ref_vcf_path(args.ref, args.population, build_num)

    if vcf is None:
        logging.info(
            "Build %s: no 1KG reference VCF for this build — "
            "VCF will be resolved after liftover to hg38.", build_num,
        )
        if args.fill_eaf:
            logging.warning(
                "--fill-eaf is not supported for build %s (positions must match "
                "the 1KG VCF build). EAF lookup will be skipped.", build_num,
            )
    else:
        if not os.path.isfile(vcf):
            logging.error("Reference VCF not found: %s", vcf)
            sys.exit(1)
        assert_bgzf(vcf, "Reference VCF")
        logging.info("Reference VCF : %s", vcf)

    # If liftover is requested, output build will be hg38.  Update build_num
    # before constructing stems/pkl_path so filenames are consistent across all
    # stages even when run_processing has not yet executed.
    if args.liftover and build_num != "38":
        build_num = "38"

    stem     = file_tag(args.gwas, args.population, input_build, build_num, added_n)
    pkl_path = os.path.join(output_loc, f"{stem}.pkl")

    # ── Stage routing ──────────────────────────────────────────────────────────
    # --only-qc is a legacy alias for --stage qc.
    stage = args.stage
    if args.only_qc:
        if stage == "all":
            stage = "qc"
            logging.info("--only-qc detected: treating as --stage qc.")
        else:
            logging.warning("--only-qc is redundant when --stage is set; ignoring --only-qc.")

    # Guard incompatible flag combinations
    _pickle_required_stages = (
        "qc", "cojo",
        "process-split",
        "process-check-ref", "process-infer-strand",
        "process-assign-rsid", "process-check-af",
        "merge",
    )
    if stage in _pickle_required_stages and args.no_pickle:
        logging.error(
            "--stage %s requires a pickle checkpoint written by a prior stage. "
            "--no-pickle is incompatible with --stage %s.", stage, stage,
        )
        sys.exit(1)

    if stage == "process-assign-rsid" and not args.dbsnp:
        logging.error(
            "--stage process-assign-rsid requires --dbsnp.  "
            "Without --dbsnp there is nothing to do in this stage."
        )
        sys.exit(1)

    # Validate --chrom range
    if args.chrom is not None:
        if not (1 <= args.chrom <= 26):
            logging.error(
                "--chrom %d is out of range. Valid values: 1–22 (autosomes), "
                "23 (X), 24 (Y), 25 (nonPAR), 26 (MT).", args.chrom,
            )
            sys.exit(1)
        _chrom_stages = {
            "process-check-ref", "process-infer-strand",
            "process-assign-rsid", "process-check-af",
        }
        if stage not in _chrom_stages:
            logging.error(
                "--chrom is only valid with per-chromosome stages: %s. "
                "Got --stage %s.", ", ".join(sorted(_chrom_stages)), stage,
            )
            sys.exit(1)

    logging.info("Pipeline stage : %s%s", stage,
                 f"  (chr {args.chrom})" if args.chrom else "")

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: preprocess
    # Load data, verify build, standardise columns.
    # --stage all  : gwas_data flows directly into the process stage below.
    # --stage preprocess : saves checkpoint and exits.
    # ══════════════════════════════════════════════════════════════════════════
    gwas_data = None
    if stage in ("all", "preprocess"):
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: preprocess")
        logging.info("=" * 60)

        logging.info("\n===== Loading GWAS summary statistics =====")
        logging.info("Loading: %s", input_path)
        _sep    = detect_separator(input_path)
        _engine = "python" if _sep == r"\s+" else "pyarrow"
        gwas_data = pd.read_csv(
            input_path, sep=_sep, header=0,
            na_values=["NA"], dtype={"CHR": "string"}, engine=_engine,
        )
        logging.info("Loaded %s variants, %d columns.",
                     f"{len(gwas_data):,}", gwas_data.shape[1])

        REFERENCE = verify_build(gwas_data, REFERENCE)

        if args.fill_eaf and vcf is not None:
            gwas_data = check_and_fill_eaf(gwas_data, vcf)
        elif args.fill_eaf:
            logging.info("EAF lookup skipped (no reference VCF for build %s).", build_num)
        else:
            logging.info("EAF lookup skipped (--fill-eaf not set).")

        if any(v is not None for v in [args.n, args.n_cases, args.n_controls]):
            gwas_data = apply_fixed_n(
                gwas_data,
                n=args.n, n_cases=args.n_cases, n_controls=args.n_controls,
                force=args.force_n,
            )

        gwas_data = correct_columns(gwas_data)
        gwas_data = standardise_columns(gwas_data)

        if args.figures:
            raw_stem = file_tag(args.gwas, args.population,
                                input_build, input_build, added_n)
            plot_raw_histograms(gwas_data, raw_stem, plots_loc)

        if stage == "preprocess":
            save_preprocess_checkpoint(
                gwas_data,
                {"reference": REFERENCE, "build_num": build_num,
                 "input_build": input_build},
                stem, output_loc,
            )
            logging.info("Stage 'preprocess' complete.")
            logging.info("Next: --stage process  (pass the same --gwas / --build / "
                         "--liftover / --output flags).")
            logging.info("=" * 70)
            logging.info("Log : %s", log_path)
            logging.info("=" * 70)
            return
        # else: gwas_data flows into the process stage below

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-normalize
    # basic_check + remove_dup + liftover.
    # --stage all / process-normalize : entry from preprocess (gwas_data in memory
    #                                   or loaded from .preprocess.parquet).
    # Saves: {stem}.normalize.pkl
    # ══════════════════════════════════════════════════════════════════════════
    gwas_obj = None
    _process_substages = {
        "process-normalize", "process-check-ref",
        "process-infer-strand", "process-assign-rsid", "process-check-af",
    }
    if stage in ("all",) or stage in _process_substages:

        if stage == "process-normalize":
            logging.info("\n" + "=" * 60)
            logging.info("STAGE: process-normalize")
            logging.info("=" * 60)
            gwas_data, meta = load_preprocess_checkpoint(stem, output_loc)
            REFERENCE   = meta.get("reference",   REFERENCE)
            build_num   = meta.get("build_num",   build_num)
            input_build = meta.get("input_build", input_build)
            stem        = file_tag(args.gwas, args.population, input_build, build_num, added_n)
            gwas_obj = make_sumstats_object(gwas_data, REFERENCE)
            del gwas_data
            gc.collect()
            gwas_obj, REFERENCE = run_normalize(
                gwas_obj, REFERENCE, args.ref, args.threads,
                args.liftover, args.population,
            )
            build_num = normalise_build(REFERENCE)
            stem      = file_tag(args.gwas, args.population, input_build, build_num, added_n)
            save_process_checkpoint(gwas_obj, stem, output_loc, "normalize")
            logging.info("Stage 'process-normalize' complete.")
            logging.info("Next: --stage process-check-ref")
            logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
            return

        elif stage == "all":
            # In --stage all, gwas_data still in memory from preprocess block
            if gwas_data is None:
                logging.error("Internal error: gwas_data not available for --stage all.")
                sys.exit(1)
            logging.info("\n" + "=" * 60)
            logging.info("STAGE: process  (all — running full run_processing)")
            logging.info("=" * 60)
            gwas_obj = make_sumstats_object(gwas_data, REFERENCE)
            del gwas_data
            gc.collect()
            gwas_obj, REFERENCE = run_processing(
                gwas_obj, REFERENCE, args.ref, vcf,
                args.threads, args.liftover, args.dbsnp,
                population=args.population,
            )
            build_num = normalise_build(REFERENCE)
            vcf       = ref_vcf_path(args.ref, args.population, build_num)
            stem      = file_tag(args.gwas, args.population, input_build, build_num, added_n)
            pkl_path  = os.path.join(output_loc, f"{stem}.pkl")
            save_raw_outputs(gwas_obj, args.gwas, args.population,
                             input_build, build_num, output_loc, added_n,
                             save_pickle=not args.no_pickle)
            if args.figures:
                plot_full_dataset(gwas_obj, args.gwas, REFERENCE,
                                  plots_loc, args.daf_max, stem)
            if args.cojo:
                write_cojo(
                    gwas_obj, args.gwas, args.population,
                    input_build, build_num, output_loc,
                    snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
                    suffix="", added_n=added_n,
                )
            # Falls through to the qc/leads/cojo blocks below

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-split
    # Split the whole-genome normalize.pkl into per-chromosome parquets.
    # Loads: {stem}.normalize.pkl
    # Saves: {stem}.chr{N}.normalize.parquet  ×  N chromosomes
    #        {stem}.chrsplit.json              (manifest)
    # Trivial resources — just DataFrame groupby + parquet writes.
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "process-split":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: process-split")
        logging.info("=" * 60)
        gwas_obj = load_process_checkpoint(stem, output_loc, "normalize")
        manifest = split_by_chrom(gwas_obj, stem, output_loc)
        logging.info("Stage 'process-split' complete. %d chromosome(s) written.",
                     len(manifest))
        logging.info("Next: --stage process-check-ref --chrom N  (SLURM array 1-26)")
        logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
        return

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-check-ref
    # check_ref + flip_allele_stats + fix_id.
    #
    # Whole-genome path (--chrom not set):
    #   Loads: {stem}.normalize.pkl
    #   Saves: {stem}.checkref.pkl
    #
    # Per-chromosome path (--chrom N):
    #   Loads: {stem}.chrN.normalize.parquet
    #   Saves: {stem}.chrN.checkref.parquet
    #   If chrN parquet missing → exits 0 (chromosome not in study).
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "process-check-ref":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: process-check-ref")
        logging.info("=" * 60)
        if args.chrom:
            df = load_chrom_parquet(stem, output_loc, args.chrom, "normalize")
            gwas_obj = make_sumstats_from_chrom_df(df, build_num)
            gwas_obj = run_check_ref(gwas_obj, build_num, args.ref)
            save_chrom_parquet(gwas_obj.data, stem, output_loc, args.chrom, "checkref")
            logging.info("Stage 'process-check-ref' (chr %d) complete.", args.chrom)
            logging.info("Next: --stage process-infer-strand --chrom %d", args.chrom)
        else:
            gwas_obj = load_process_checkpoint(stem, output_loc, "normalize")
            gwas_obj = run_check_ref(gwas_obj, build_num, args.ref)
            save_process_checkpoint(gwas_obj, stem, output_loc, "checkref")
            logging.info("Stage 'process-check-ref' complete.")
            logging.info("Next: --stage process-infer-strand")
        logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
        return

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-infer-strand
    # infer_strand2 + flip_allele_stats.
    #
    # Whole-genome path (--chrom not set):
    #   Loads: {stem}.checkref.pkl
    #   Saves: {stem}.inferstrand.pkl
    #   High memory — full 1KG VCF sweep.
    #
    # Per-chromosome path (--chrom N):
    #   Loads: {stem}.chrN.checkref.parquet
    #   Saves: {stem}.chrN.inferstrand.parquet
    #   Memory proportional to chromosome variant count (~1/22 of full sweep).
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "process-infer-strand":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: process-infer-strand")
        logging.info("=" * 60)
        _vcf = ref_vcf_path(args.ref, args.population, build_num)
        if args.chrom:
            df = load_chrom_parquet(stem, output_loc, args.chrom, "checkref")
            gwas_obj = make_sumstats_from_chrom_df(df, build_num)
            gwas_obj = run_infer_strand(gwas_obj, args.ref, _vcf, args.threads)
            save_chrom_parquet(gwas_obj.data, stem, output_loc, args.chrom, "inferstrand")
            _next = "process-assign-rsid" if args.dbsnp else "process-check-af"
            logging.info("Stage 'process-infer-strand' (chr %d) complete.", args.chrom)
            logging.info("Next: --stage %s --chrom %d", _next, args.chrom)
        else:
            gwas_obj = load_process_checkpoint(stem, output_loc, "checkref")
            gwas_obj = run_infer_strand(gwas_obj, args.ref, _vcf, args.threads)
            save_process_checkpoint(gwas_obj, stem, output_loc, "inferstrand")
            _next = "process-assign-rsid" if args.dbsnp else "process-check-af"
            logging.info("Stage 'process-infer-strand' complete.")
            logging.info("Next: --stage %s", _next)
        logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
        return

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-assign-rsid
    # assign_rsid via dbSNP VCF sweep.
    # Only reached when --dbsnp is set (guarded above).
    #
    # Whole-genome path (--chrom not set):
    #   Loads: {stem}.inferstrand.pkl
    #   Saves: {stem}.assignrsid.pkl
    #   Extreme memory — full dbSNP VCF sweep.
    #
    # Per-chromosome path (--chrom N):
    #   Loads: {stem}.chrN.inferstrand.parquet
    #   Saves: {stem}.chrN.assignrsid.parquet
    #   Memory proportional to chromosome variant count.
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "process-assign-rsid":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: process-assign-rsid")
        logging.info("=" * 60)
        if args.chrom:
            df = load_chrom_parquet(stem, output_loc, args.chrom, "inferstrand")
            gwas_obj = make_sumstats_from_chrom_df(df, build_num)
            gwas_obj = run_assign_rsid(gwas_obj, build_num, args.ref, args.threads)
            save_chrom_parquet(gwas_obj.data, stem, output_loc, args.chrom, "assignrsid")
            logging.info("Stage 'process-assign-rsid' (chr %d) complete.", args.chrom)
            logging.info("Next: --stage process-check-af --chrom %d", args.chrom)
        else:
            gwas_obj = load_process_checkpoint(stem, output_loc, "inferstrand")
            gwas_obj = run_assign_rsid(gwas_obj, build_num, args.ref, args.threads)
            save_process_checkpoint(gwas_obj, stem, output_loc, "assignrsid")
            logging.info("Stage 'process-assign-rsid' complete.")
            logging.info("Next: --stage process-check-af")
        logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
        return

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: process-check-af
    # check_af2 — frequency concordance against 1KG VCF.
    #
    # Whole-genome path (--chrom not set):
    #   Loads: {stem}.assignrsid.pkl (--dbsnp) or {stem}.inferstrand.pkl
    #   Saves: final raw outputs — {stem}.pkl + .parquet + .tsv.gz
    #   High memory — full 1KG VCF sweep.
    #   Next stage: qc (or merge if using the chr-split path).
    #
    # Per-chromosome path (--chrom N):
    #   Loads: {stem}.chrN.assignrsid.parquet (--dbsnp) or {stem}.chrN.inferstrand.parquet
    #   Saves: {stem}.chrN.checkaf.parquet
    #   Terminal per-chr stage — merge follows after all chromosomes complete.
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "process-check-af":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: process-check-af")
        logging.info("=" * 60)
        _vcf = ref_vcf_path(args.ref, args.population, build_num)
        _prev_suffix = "assignrsid" if args.dbsnp else "inferstrand"
        if args.chrom:
            df = load_chrom_parquet(stem, output_loc, args.chrom, _prev_suffix)
            gwas_obj = make_sumstats_from_chrom_df(df, build_num)
            gwas_obj = run_check_af(gwas_obj, _vcf, args.threads)
            save_chrom_parquet(gwas_obj.data, stem, output_loc, args.chrom, "checkaf")
            logging.info("Stage 'process-check-af' (chr %d) complete.", args.chrom)
            logging.info("After all chromosomes: --stage merge")
        else:
            gwas_obj = load_process_checkpoint(stem, output_loc, _prev_suffix)
            gwas_obj = run_check_af(gwas_obj, _vcf, args.threads)
            pkl_path = os.path.join(output_loc, f"{stem}.pkl")
            # Whole-genome terminal process sub-stage: write final raw outputs
            save_raw_outputs(gwas_obj, args.gwas, args.population,
                             input_build, build_num, output_loc, added_n,
                             save_pickle=not args.no_pickle)
            if args.figures:
                plot_full_dataset(gwas_obj, args.gwas, build_num,
                                  plots_loc, args.daf_max, stem)
            if args.cojo:
                write_cojo(
                    gwas_obj, args.gwas, args.population,
                    input_build, build_num, output_loc,
                    snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
                    suffix="", added_n=added_n,
                )
            logging.info("Stage 'process-check-af' complete.")
            logging.info("Next: --stage qc")
        logging.info("=" * 70); logging.info("Log : %s", log_path); logging.info("=" * 70)
        return

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: merge
    # Concatenate all per-chromosome checkaf parquets, then run QC + plots +
    # leads + COJO.  Replaces the qc + cojo stages in the chr-split path.
    # Loads: {stem}.chrsplit.json  +  {stem}.chr{N}.checkaf.parquet × N
    # Saves: same final outputs as process-check-af (raw) + qc + cojo stages
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "merge":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: merge")
        logging.info("=" * 60)
        run_merge(stem, output_loc, REFERENCE, build_num, input_build, args)
        logging.info("Stage 'merge' complete.")
        logging.info("=" * 70)
        logging.info("Output : %s", output_loc)
        logging.info("Log    : %s", log_path)
        logging.info("=" * 70)
        return

    # Ensure stem is always defined for both the normal path and the qc/cojo
    # stages that enter below without running the process block.
    stem = file_tag(args.gwas, args.population, input_build, build_num, added_n)

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: qc
    # Apply QC filters, save QC outputs, generate QC plots, extract leads.
    # --stage all : gwas_obj is still in memory from the process stage above.
    # --stage qc  : loads the raw pickle written by --stage process.
    # ══════════════════════════════════════════════════════════════════════════
    gwas_obj_qc = None
    run_qc = (args.qc or stage == "qc") and stage != "cojo"
    if run_qc:
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: qc")
        logging.info("=" * 60)

        if stage == "qc":
            logging.info("[qc] Loading pickle: %s", pkl_path)
            gwas_obj = gl.load_pickle(pkl_path)
            if gwas_obj is None:
                logging.error(
                    "[qc] Failed to load pickle — file not found or unreadable:\n"
                    "  %s\n"
                    "Ensure --stage process has been run with the same "
                    "--gwas / --build / --liftover / --output flags.",
                    pkl_path,
                )
                sys.exit(1)
            # Optionally regenerate full-dataset plots from the loaded pickle.
            if args.figures:
                plot_full_dataset(gwas_obj, args.gwas, build_num,
                                  plots_loc, args.daf_max, stem)

        gwas_obj_qc = apply_qc(
            gwas_obj,
            args.eaf_min, args.beta_max, args.se_max,
            args.info_min, args.daf_max,
        )
        # Free the unfiltered object — gwas_obj_qc is now the working copy.
        del gwas_obj
        gc.collect()

        save_qc_outputs(gwas_obj_qc, args.gwas, args.population,
                        input_build, build_num, output_loc, added_n,
                        save_pickle=not args.no_pickle)

        if args.figures:
            plot_qc_dataset(gwas_obj_qc, args.gwas, build_num,
                            plots_loc, args.daf_max, stem)

        if args.cojo:
            write_cojo(
                gwas_obj_qc, args.gwas, args.population,
                input_build, build_num, output_loc,
                snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
                suffix="qc", added_n=added_n,
            )

        if args.leads:
            extract_leads(gwas_obj_qc, args.gwas, args.population,
                          input_build, build_num, output_loc, True, added_n)

        if stage == "qc":
            logging.info("Stage 'qc' complete.")
            logging.info("=" * 70)
            logging.info("Output : %s", output_loc)
            logging.info("Log    : %s", log_path)
            logging.info("=" * 70)
            return

    # ── Leads (--stage all, QC disabled) ──────────────────────────────────────
    # When QC ran, leads were already extracted inside the qc block above.
    # When QC did not run, extract from gwas_obj which is still in scope.
    if stage == "all" and args.leads and not run_qc:
        if gwas_obj is not None:
            extract_leads(gwas_obj, args.gwas, args.population,
                          input_build, build_num, output_loc, False, added_n)
        else:
            logging.warning("--leads requested but gwas_obj is not available; skipping.")

    # ══════════════════════════════════════════════════════════════════════════
    # STAGE: cojo  (standalone)
    # Write a COJO file from an existing pickle.  Prefers the QC pickle if
    # present; falls back to the raw pickle.
    # ══════════════════════════════════════════════════════════════════════════
    if stage == "cojo":
        logging.info("\n" + "=" * 60)
        logging.info("STAGE: cojo")
        logging.info("=" * 60)

        qc_pkl_path = os.path.join(output_loc, f"{stem}.qc.pkl")
        if os.path.isfile(qc_pkl_path):
            logging.info("[cojo] Loading QC pickle: %s", qc_pkl_path)
            cojo_obj    = gl.load_pickle(qc_pkl_path)
            cojo_suffix = "qc"
        elif os.path.isfile(pkl_path):
            logging.info("[cojo] Loading raw pickle: %s", pkl_path)
            cojo_obj    = gl.load_pickle(pkl_path)
            cojo_suffix = ""
        else:
            logging.error(
                "[cojo] No pickle found at:\n  %s\n  %s\n"
                "Run --stage process (and optionally --stage qc) first.",
                pkl_path, qc_pkl_path,
            )
            sys.exit(1)

        write_cojo(
            cojo_obj, args.gwas, args.population,
            input_build, build_num, output_loc,
            snpid_fmt=args.cojo_id, add_pos=args.cojo_pos,
            suffix=cojo_suffix, added_n=added_n,
        )
        logging.info("Stage 'cojo' complete.")
        logging.info("=" * 70)
        logging.info("Output : %s", output_loc)
        logging.info("Log    : %s", log_path)
        logging.info("=" * 70)
        return

    # ── Done ───────────────────────────────────────────────────────────────────
    logging.info("=" * 70)
    logging.info("Pipeline complete.")
    logging.info("  Output : %s", output_loc)
    logging.info("  Plots  : %s", plots_loc)
    logging.info("  Log    : %s", log_path)
    logging.info("=" * 70)

if __name__ == "__main__":
    main()
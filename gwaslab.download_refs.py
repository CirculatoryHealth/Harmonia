#!/usr/bin/env python3
"""
gwaslab.download_refs.py — Download gwaslab reference files.

Usage
-----
    python gwaslab.download_refs.py [--ref-dir DIR]

Default ref dir: /hpc/dhl_ec/data/references/gwaslab/

What this downloads
-------------------
All 1KG population VCFs (+ auto-downloaded .tbi index) for AFR, EAS, AMR, SAS
on both hg19 and hg38.

All other reference files are already present at the default ref dir and are
listed here as comments for completeness. Uncomment any entry to (re-)download it.

Note: recombination maps and GTF files are downloaded automatically by gwaslab
at runtime when needed; pre-downloading them here is optional but useful on HPC
nodes without outbound internet during jobs.

Reference
---------
https://cloufield.github.io/gwaslab/Download/
"""

import argparse
import gwaslab as gl
from gwaslab.bd.bd_download import (
    set_default_directory,
    get_default_directory,
    download_ref,
    check_downloaded_ref,
    get_path,
)

DEFAULT_REF_DIR = "/hpc/dhl_ec/data/references/gwaslab/"

# ---------------------------------------------------------------------------
# Reference file keywords.
# Uncomment any line to (re-)download that file.
# .tbi index files are fetched automatically alongside VCFs.
# ---------------------------------------------------------------------------
# fmt: off
TO_DOWNLOAD = [

    # ── 1KG population VCFs — hg19 (1KGp3v5 low-coverage) ──────────────────
    # "1kg_eur_hg19",             # already present: EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
    # "1kg_pan_hg19",             # already present: PAN.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
    "1kg_afr_hg19",
    "1kg_eas_hg19",
    "1kg_amr_hg19",
    "1kg_sas_hg19",

    # ── 1KG population VCFs — hg38 (1KGp 30x) ──────────────────────────────
    # "1kg_eur_hg38",             # already present: EUR.ALL.split_norm_af.1kg_30x.hg38.vcf.gz
    # "1kg_pan_hg38",             # already present: PAN.ALL.split_norm_af.1kg_30x.hg38.vcf.gz
    "1kg_afr_hg38",
    "1kg_eas_hg38",
    "1kg_amr_hg38",
    "1kg_sas_hg38",

    # ── HapMap3 EAF (pan-ancestry allele frequencies) ───────────────────────
    # "1kg_hm3_hg19_eaf",         # already present: PAN.hapmap3.hg19.EAF.tsv.gz
    # "1kg_hm3_hg38_eaf",         # already present: PAN.hapmap3.hg38.EAF.tsv.gz

    # ── 1KG SNPID→rsID conversion tables (no VCF sweep needed) ─────────────
    # "1kg_dbsnp151_hg19_auto",   # already present: 1kg_dbsnp151_hg19_auto.txt.gz
    # "1kg_dbsnp151_hg38_auto",   # already present: 1kg_dbsnp151_hg38_auto.txt.gz

    # ── dbSNP VCFs (!very large, NCBI FTP) ──────────────────────────────────
    # "dbsnp_v151_hg19",          # already present: 00-All.vcf.gz (hg19)
    # "dbsnp_v151_hg38",          # v151 hg38
    # "dbsnp_v157_hg19",          # already present: GCF_000001405.25.gz
    # "dbsnp_v157_hg38",          # already present: GCF_000001405.40.gz

    # ── Reference genome FASTA (UCSC, !large) ───────────────────────────────
    # "ucsc_genome_hg19",         # already present: hg19.fa.gz
    # "ucsc_genome_hg38",         # already present: hg38.fa.gz

    # ── Recombination maps (auto-downloaded by gwaslab at runtime) ───────────
    # "recombination_hg19",       # already present: recombination_hg19.tar.gz
    # "recombination_hg38",       # already present: recombination_hg38.tar.gz

    # ── Gene annotation GTFs (auto-downloaded by gwaslab at runtime) ────────
    # "ensembl_hg19_gtf",         # already present: Homo_sapiens.GRCh37.87.chr.gtf.gz
    # "ensembl_hg38_gtf",         # already present: Homo_sapiens.GRCh38.109.chr.gtf.gz
    # "refseq_hg19_gtf",          # already present: GRCh37_latest_genomic.gtf.gz
    # "refseq_hg38_gtf",          # already present: GRCh38_latest_genomic.gtf.gz
]
# fmt: on


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ref-dir", default=DEFAULT_REF_DIR,
                        help=f"Reference directory (default: {DEFAULT_REF_DIR})")
    args = parser.parse_args()

    print(f"Reference directory : {args.ref_dir}")
    set_default_directory(args.ref_dir)

    for keyword in TO_DOWNLOAD:
        print(f"\n{'─' * 60}")
        print(f"Downloading: {keyword}")
        download_ref(keyword)
        path = get_path(keyword, verbose=False)
        print(f"  → {path}")

    print(f"\n{'═' * 60}")
    print("Done. Downloaded references:")
    check_downloaded_ref()


if __name__ == "__main__":
    main()

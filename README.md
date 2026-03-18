[gwas2cojo](https://github.com/CirculatoryHealth/gwas2cojo)<img align="right" height="200" src=logo/fulllogo_transparent.png>
============
[![Languages](https://skillicons.dev/icons?i=bash,py)](https://skillicons.dev) 


# 📑 Introduction
This suite of scripts offers two standalone pipelines to process genome-wide association study (GWAS) summary statistics datasets relative to a reference and parse them into a standardised format for downstream analyses. 

* The first script, `gwas2cojo.py`, is a public python script that aligns a public GWAS dataset to a genetic reference, to enable large scale cross dataset comparisons with public available GWAS datasets. This relies on the 1000G phase 3 reference dataset for Europeans, but can be used with any reference dataset with the appropriate columns. It generates a [COJO]-compatible file that can be used for many post-GWAS analyses, including SMR. It is the recommended tool for legacy datasets using b37 and hundreds of thousands of variants (< ±7 million).

* The second script, `gwaslab.process.py`, is a more recent and comprehensive pipeline built on the [GWASLab](https://github.com/Cloufield/gwaslab) library. It relies on the human reference, dbSNP, and 1000G phase 3 data. It is the recommended successor to `gwas2cojo.py` for new datasets using b38 and tens of million variants.

## ⚙️ Requirements

* Python 3.12
* Tested on Linux _Rocky Linux release 8.10 (Green Obsidian)_ and macOS _Tahoe (26.4 Beta (25E5233c))_
* `bcftools` ≥1.17 (used by GWASLab VCF sweep steps; installed via bioconda)
* See `environment.yml` for the full dependency list

## 🧩 Installation

We recommend using `mamba` (a faster drop-in replacement for `conda`).

### Option A — from `environment.yml` *(recommended)*

```bash
mamba env create -f environment.yml
mamba activate gwas2cojo
```

This installs Python 3.12, `bcftools` (bioconda), and all pip dependencies in one step.

### Option B — manual create + pip install

If you prefer to build the environment by hand:

```bash
mamba create --name gwas2cojo python=3.12
mamba activate gwas2cojo
pip install \
    "numpy>=1.21.2,<2" pyliftover tqdm "adjusttext==0.8" \
    "matplotlib>=3.8,<3.9" "pandas>=1.3,!=1.5" "pysam==0.22.1" \
    "scikit-allel>=1.3.5" "scipy>=1.12" "seaborn>=0.12" \
    "h5py>=3.10.0" pyarrow "polars>=1.27.0" \
    "sumstats-liftover==1.1.0" "jupyter==1.0.0" \
    gwaslab bcftools
```

> **Note:** `bcftools` is also available as a conda package from bioconda and is preferred over the pip wrapper on most HPC systems:
> ```bash
> mamba install -c bioconda bcftools
> ```

### 2️⃣ Verify installation

Check that the core modules load correctly:

```bash
python -c "import numpy, gwaslab, pyliftover, polars; print('OK')"
```

If you see `OK`, the environment is ready.

### 3️⃣ (Optional) Test the scripts

Confirm that `gwas2cojo.py` runs:

```bash
python gwas2cojo.py --help
```

Confirm that `gwaslab.process.py` runs:

```bash
python gwaslab.process.py --help
```

For a minimal parsing test without running the full pipeline:

```bash
python gwas2cojo.py --header-only --gwas example_gwas.txt.gz
```

## 🔧 Troubleshooting

### 🧱 Dependency conflicts

If `pip install` reports conflicts, try installing inside the activated conda environment and letting conda resolve system libraries first:

```bash
mamba install -c conda-forge -c bioconda bcftools pysam h5py
pip install "numpy>=1.21.2,<2" gwaslab "polars>=1.27.0" ...
```

### 🐍 Prefer `conda` instead of `mamba`?

```bash
conda env create -f environment.yml
conda activate gwas2cojo
```

---

# 🔭 gwas2cojo.py — Aligns GWAS to cojo-format

`gwas2cojo.py` is a public python script that aligns a public GWAS dataset to a genetic reference, to enable large scale cross dataset comparisons with public available GWAS datasets. Among others, it tries to deal with different dataformats and different genome builds.
Most importantly, it aligns the variant notation such that swapped, translated, wrong or ambiguous ambivalent alleles are corrected or removed.

## 🧩 Output format

It generates a [COJO]-compatible file:

```
SNP       A1  A2  freq    b       se      p       n
rs1001    A   G   0.8493  0.0024  0.0055  0.6653  129850
rs1002    C   G   0.03606 0.0034  0.0115  0.7659  129799
rs1003    A   C   0.5128  0.045   0.038   0.2319  129830
```

[COJO]: https://yanglab.westlake.edu.cn/software/smr/#Overview

## 🚀 Usage

The usage is intuitive and many options are provided.

```
python3 gwas2cojo.py \
        --gen:build     hg19 \
        --gen           1kGp3.ref.1maf.nonbia.sumstats.gz \
        --gwas          "${gwasfilename}" \
        --report        "out/${name}.report" \
        --gen:ident     ID \
        --gen:chr       CHROM \
        --gen:other     REF \
        --gen:effect    ALT \
        --gen:eaf       AF \
        --out           "out/${name}.cojo"
```

## 📜 List of options

A full list of options is given below.

```
$ python3 ./gwas2cojo.py -h
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        CONVERT GWAS FOR SMR


* Written by         : Lennart Landsmeer | lennart[at]landsmeer[dot]email
* Suggested for by   : Sander W. van der Laan | s.w.vanderlaan[at]gmail[dot]com
* Last update        : 2024-04-20
* Name               : gwas2cojo
* Version            : v1.4.3

* Description        : Converts a given set of summary statistics from genome-wide association studies
                       (GWAS) to the GWAS-COJO format used by Summarized-data Mendelian Randomization
                       (SMR). This format is also usable for many other post-GWAS analyses.
                       A reference, e.g. 1000G phase 3, is used to map GWAS SumStats to.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Start: 2025-11-07 04:53:50.675318
usage: gwas2cojo.py [-h] [-o cojo] [-r txt] [-rr] [-g file.stats.gz] --gwas file.txt.gz. [--header-only] [--output-pos] [--fmid MID] [--fclose CLOSE]
                    [--ignore-indels] [--gwas:effect COLUMN] [--gwas:other COLUMN] [--gwas:freq COLUMN] [--gwas:beta COLUMN] [--gwas:or COLUMN]
                    [--gwas:se COLUMN] [--gwas:p COLUMN] [--gwas:chr-bp COLUMN] [--gwas:chr COLUMN] [--gwas:bp COLUMN] [--gwas:build BUILDID] [--gwas:n COLUMNS]
                    [--gwas:sep DELIMITER] [--gwas:header:remove STRING] [--gwas:default:p VALUE] [--gwas:default:beta VALUE] [--gwas:default:se VALUE]
                    [--gwas:default:chr VALUE] [--gwas:default:n VALUE] [--gen:ident COLUMN] [--gen:chr COLUMN] [--gen:bp COLUMN] [--gen:effect COLUMN]
                    [--gen:other COLUMN] [--gen:eaf COLUMN] [--gen:oaf COLUMN] [--gen:maf COLUMN] [--gen:minor COLUMN] [--gen:build COLUMN]

optional arguments:
  -h, --help            show this help message and exit
  -o cojo, --out cojo   Output .cojo file.
  -r txt, --report txt  Report discarded variants here.
  -rr, --report-ok      Report all decisions made. Warning: very verbose
  -g file.stats.gz, --gen file.stats.gz
                        Genetic reference data. Could be an in-house GWAS or a reference dataset (e.g. 1000G phase1, phase3, etc.) with the following columns:
                        CHROM, POS, ID, REF, ALT, CHROM:POS:REF:ALT, AF, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF.
  --gwas file.txt.gz.   GWAS summary statistics location.
  --header-only         Exit after reading GWAS header. Useful for testing whether a file is readable by this program.
  --output-pos          Write chr and bp columns to output file (breaks COJO format assumptions)

filter snps:
  --fmid MID            Ambivalent variants are ambiguous when effect frequency is between 0.5-MID and 0.5+MID. Set to 0 to prevent discarding. Default is 0.05.
  --fclose CLOSE        Fequencies are considered close when their difference is less than CLOSE. Set to 1 to prevent discarding. Default is 0.1
  --ignore-indels       Should insertions and deletions be ignored? Only SNPs are retained.

gwas header:
  --gwas:effect COLUMN  Effect/Risk/Coded/Minor allele column name.
  --gwas:other COLUMN   Non-effect/Other/Major allele column name.
  --gwas:freq COLUMN    Effect/Risk/Coded/Minor allele frequency column name.
  --gwas:beta COLUMN    Log-odds column name [beta/effect], relative to effect/risk/coded/minor allele.
  --gwas:or COLUMN      Odds-ratio column name [OR], relative to effect/risk/coded/minor allele.
  --gwas:se COLUMN      Log-odds standard error column name.
  --gwas:p COLUMN       P-value column name.
  --gwas:chr-bp COLUMN  Position column name when encoded as chr:pos.
  --gwas:chr COLUMN     Chromosome column name.
  --gwas:bp COLUMN      Chromosomal position column name.
  --gwas:build BUILDID  hg18 or b36, hg19 or b37, etc.
  --gwas:n COLUMN(S)    Column name(s) of the sample counts. Separated by commas. If multiple colums are specified, their sum is stored.
  --gwas:sep DELIMITER  Delimiter character. Defaults to any whitespace character
  --gwas:header:remove STRING
                        Remove this string from header before processing

gwas default values:
  --gwas:default:p VALUE
  --gwas:default:beta VALUE
  --gwas:default:se VALUE
  --gwas:default:chr VALUE
  --gwas:default:n VALUE

genetic data header:
  --gen:ident COLUMN    Column name of variant identifier (e.g. rsid).
  --gen:chr COLUMN      Column name of chromosome.
  --gen:bp COLUMN       Column name of bp position.
  --gen:effect COLUMN   Column name of effect allele.
  --gen:other COLUMN    Column name of non-effect allele.
  --gen:eaf COLUMN      Column name of effect allele frequency.
  --gen:oaf COLUMN      Column name of non-effect allele frequency.
  --gen:maf COLUMN      Column name of minor allele frequency.
  --gen:minor COLUMN    Column name of minor allele. When used in combination with maf, it is used to find the effect allele frequency.
  --gen:build COLUMN    Genetic reference data build. Defaults to hg19
```

See also [this link] for more background and an additional explanation.

[this link]: https://blog.llandsmeer.com/tech/2019/12/28/gwas2cojo.html

## 🧠 Reference datasets

You will need a reference to map the data to. You can create your own, or use the one we created [one based](https://blog.llandsmeer.com/1kGp3.ref.1maf.nonbia.sumstats.gz) on the 1000G phase 3 data for Europeans. This is filtered based on MAF>1% and excludes non-bi-allelic and duplicate variants.


# 🔬 gwaslab.process.py — GWASLab Processing Pipeline

`gwaslab.process.py` is a standalone pipeline built on the [GWASLab](https://github.com/Cloufield/gwaslab) library. It takes raw GWAS summary statistics and runs them through a fully automated processing chain: standardisation, strand inference, build liftover, dbSNP annotation, allele-frequency validation, QC filtering, and output generation in multiple formats (`pickle`, `parquet`, `tsv.gz`, `cojo.gz`). It is the recommended successor to `gwas2cojo.py` for new datasets using b38 and tens of million variants.


## 🗺️ Pipeline overview

Steps run in order; individual steps can be toggled with the flags described below. Each stage can be run as an independent SLURM job using `--stage`; see [HPC: Staged SLURM submission](#-hpc-staged-slurm-submission) below.

| Step | Description | Stage (`--stage`) | Toggle |
|------|-------------|-------------------|--------|
| 1 | Load and parse input file; auto-detect column names | `preprocess` | always |
| 2 | Verify / normalise genome build | `preprocess` | always |
| 3 | Correct P-values, SE, and N columns | `preprocess` | always |
| 4 | Standardise to GWASLab `Sumstats` object | `preprocess` | always |
| 5 | Plot raw input histograms | `preprocess` | `--figures` |
| 6 | `basic_check` — flag malformed variants | `process-normalize` | always |
| 7 | `remove_dup` — drop duplicate variants | `process-normalize` | always |
| 8 | Liftover (hg18→hg38 or hg19→hg38) | `process-normalize` | `--liftover` |
| 9 | `check_ref` against hg38 FASTA | `process-check-ref` | always |
| 10 | `flip_allele_stats` — correct BETA/EAF for reference-flipped alleles | `process-check-ref` | always |
| 11 | `fix_id` — normalise variant IDs | `process-check-ref` | always |
| 12 | `infer_strand2` — resolve strand ambiguity from 1KG VCF (full sweep) | `process-infer-strand` | always |
| 13 | `flip_allele_stats` — correct stats for strand-resolved variants | `process-infer-strand` | always |
| 14 | `assign_rsid` — annotate with dbSNP rsIDs (full dbSNP VCF sweep) | `process-assign-rsid` | `--dbsnp` |
| 15 | `check_af2` — allele-frequency concordance vs 1KG VCF | `process-check-af` | always |
| 16 | Save raw outputs (pickle + parquet + TSV.GZ) | `process-check-af` | always |
| 17 | Manhattan, QQ, and DAF plots (pre-QC) | `process-check-af` | `--figures` |
| 18 | QC filter (EAF, BETA, SE, INFO, HWE, MAC, DAF) | `qc` | `--qc` |
| 19 | Save QC outputs | `qc` | `--qc` |
| 20 | Manhattan, QQ, and DAF plots (post-QC) | `qc` | `--figures` + `--qc` |
| 21 | Extract genome-wide significant lead SNPs | `qc` | `--leads` |
| 22 | Write COJO-format file(s) | `cojo` | `--cojo` |

## 🚀 Usage

### Minimal — full default pipeline (hg19 input)

```bash
python3 gwaslab.process.py \
    --gwas    MyStudy \
    --input   MyStudy.parsed.txt.gz \
    --directory /data/gwas/MyStudy \
    --ref     /data/references/gwaslab \
    --output  /data/results/MyStudy \
    --population EUR \
    --build   19 \
    --liftover --dbsnp --qc --figures --leads \
    --cojo --cojo-pos --cojo-id rsid
```

### hg18 input with all options

```bash
python3 gwaslab.process.py \
    --gwas    CAD_SCHUNKERT \
    --input   cardiogram_gwas_results_edited.txt.gz \
    --directory /data/gwas/CARDIoGRAM \
    --ref     /data/references/gwaslab \
    --output  /data/results/CAD_SCHUNKERT \
    --population EUR \
    --build   18 \
    --liftover --dbsnp --qc --figures --leads --fill-eaf \
    --cojo --cojo-pos --cojo-id rsid \
    --threads 4
```

### Run a single stage (staged mode)

Use `--stage` to run only one part of the pipeline. Pass the same `--gwas`, `--build`, `--liftover`, and `--output` flags to every stage so file stems match and checkpoints can be found.

```bash
# Stage 1 — light: load + standardise (saves .preprocess.parquet)
python3 gwaslab.process.py --gwas CAD_SCHUNKERT --input cardiogram_gwas_results_edited.txt.gz \
    --directory /data/gwas/CARDIoGRAM --ref /data/references/gwaslab \
    --output /data/results/CAD_SCHUNKERT --population EUR --build 18 \
    --liftover --dbsnp --qc --figures --leads --cojo --cojo-pos --cojo-id rsid \
    --stage preprocess

# Stage 2 — light: basic checks + liftover (saves .normalize.pkl)
python3 gwaslab.process.py ... --stage process-normalize

# Stage 3 — medium: reference check + flip (saves .checkref.pkl)
python3 gwaslab.process.py ... --stage process-check-ref

# Stage 4 — heavy: 1KG strand inference (saves .inferstrand.pkl)
python3 gwaslab.process.py ... --stage process-infer-strand

# Stage 5 — heaviest: dbSNP rsID sweep (saves .assignrsid.pkl)
python3 gwaslab.process.py ... --stage process-assign-rsid

# Stage 6 — heavy: AF check + save final raw outputs (.pkl / .parquet / .tsv.gz)
python3 gwaslab.process.py ... --stage process-check-af

# Stage 7 — medium: QC filter + plots + leads
python3 gwaslab.process.py ... --stage qc

# Stage 8 — light: write COJO file
python3 gwaslab.process.py ... --stage cojo
```

> **Important:** pass the same `--gwas`, `--population`, `--build`, `--liftover`, and `--output` flags to every stage so file stems match.

### Resume QC from a saved pickle

If the pipeline completed the processing stages but failed later (e.g. during plotting), resume from the pickle without repeating the expensive VCF sweeps:

```bash
python3 gwaslab.process.py \
    --gwas    CAD_SCHUNKERT \
    --input   cardiogram_gwas_results_edited.txt.gz \
    --directory /data/gwas/CARDIoGRAM \
    --ref     /data/references/gwaslab \
    --output  /data/results/CAD_SCHUNKERT \
    --population EUR --build 18 --liftover \
    --qc --figures --leads \
    --stage qc          # replaces the old --only-qc flag
```

## 📜 Full argument reference

### Required

| Argument | Description |
|----------|-------------|
| `--gwas NAME` | Study / phenotype name; used in all output filenames |
| `--input FILE` | Input summary-statistics file (TSV or TSV.GZ) |
| `--ref DIR` | Reference-files directory (FASTA, 1KG VCFs, dbSNP VCFs) |

### Path arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--directory DIR` | `.` | Directory containing the input file |
| `--output DIR` | `<directory>/<gwas>/GWASCatalog` | Output directory |

### Population / build

| Argument | Default | Description |
|----------|---------|-------------|
| `--population` | `EUR` | Ancestry code: `EUR` `AFR` `EAS` `AMR` `SAS` `PAN` |
| `--build` | `19` | Input genome build: `18` / `19` / `38` (and aliases `hg18`, `GRCh37`, etc.) |

### Pipeline step toggles

| Argument | Description |
|----------|-------------|
| `--liftover` | Lift coordinates to hg38 (hg18→hg38 or hg19→hg38). Requires `hg18ToHg38.over.chain.gz` in `--ref` for build 18. |
| `--dbsnp` | Assign rsIDs from a dbSNP VCF |
| `--qc` | Apply QC filters and save a filtered output |
| `--only-qc` | *(Deprecated — use `--stage qc`)* Reload from an existing pickle and run QC / plots / leads only |
| `--fill-eaf` | Look up missing EAF values from the 1KG reference VCF |
| `--figures` | Generate diagnostic plots (Manhattan, QQ, DAF, histograms) |
| `--leads` | Extract genome-wide significant lead SNPs (p < 5×10⁻⁸) |
| `--no-pickle` | Skip writing `.pkl` files (reduces peak memory and disk; incompatible with `--stage qc` and `--stage cojo`) |

### Stage control

| Argument | Description |
|----------|-------------|
| `--stage all` | Run the full pipeline end-to-end *(default)* |
| `--stage preprocess` | Load + standardise → `.preprocess.parquet` + `.preprocess.json` |
| `--stage process-normalize` | `basic_check` + `remove_dup` + liftover → `.normalize.pkl` |
| `--stage process-split` | Split `normalize.pkl` into per-chromosome parquets → `chr{N}.normalize.parquet` × N + `chrsplit.json` *(trivial)* |
| `--stage process-check-ref` | `check_ref` + `flip_allele_stats` + `fix_id`; whole-genome → `.checkref.pkl` or per-chr with `--chrom N` → `.chr{N}.checkref.parquet` |
| `--stage process-infer-strand` | `infer_strand2` + `flip_allele_stats`; whole-genome → `.inferstrand.pkl` or per-chr with `--chrom N` |
| `--stage process-assign-rsid` | `assign_rsid` via dbSNP VCF; whole-genome → `.assignrsid.pkl` or per-chr with `--chrom N` *(requires `--dbsnp`)* |
| `--stage process-check-af` | `check_af2`; whole-genome → final raw outputs `.pkl/.parquet/.tsv.gz` or per-chr with `--chrom N` → `.chr{N}.checkaf.parquet` |
| `--stage merge` | Concat all `chr{N}.checkaf.parquet` → raw outputs + QC + plots + leads + COJO *(per-chr path)* |
| `--stage qc` | QC filter → `.qc.pkl/.parquet/.tsv.gz`; loads prior pickle *(whole-genome path only)* |
| `--stage cojo` | Write COJO file from existing pickle *(whole-genome path only)* |
| `--chrom N` | Chromosome task ID 1–26 (23=X, 24=Y, 25=nonPAR, 26=MT). Set automatically by SLURM array jobs. If the chromosome is absent the stage exits 0 gracefully. |

### COJO output

| Argument | Default | Description |
|----------|---------|-------------|
| `--cojo` | off | Write a COJO-format file (`SNP A1 A2 freq b se p n`) |
| `--cojo-id` | `chrpos` | SNP column format: `chrpos` (CHR:BP) or `rsid` |
| `--cojo-pos` | off | Add CHR and BP columns after the SNP column |

### Fixed sample sizes

Useful when N is absent from the input file.

| Argument | Description |
|----------|-------------|
| `--n INT` | Total sample size |
| `--n-cases INT` | Number of cases |
| `--n-controls INT` | Number of controls |
| `--force-n` | Overwrite existing N columns with the supplied values |

### QC thresholds

All thresholds apply only when `--qc` is enabled.

| Argument | Default | Description |
|----------|---------|-------------|
| `--eaf-min F` | `0.005` | EAF lower bound (upper = 1 − F) |
| `--beta-max F` | `5.0` | Maximum \|BETA\| |
| `--se-max F` | `5.0` | Maximum SE |
| `--info-min F` | `0.4` | Minimum INFO / imputation quality score |
| `--hwe-min F` | `1e-3` | Minimum HWE p-value |
| `--mac-min N` | `30` | Minimum minor allele count |
| `--daf-max F` | `0.12` | Maximum \|DAF\| vs 1KG reference; set to `0` to disable |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads N` | `4` | Threads for parallelised gwaslab steps |

## 🗂️ Output files

All output files share a common stem:

```
{GWAS}.{POPULATION}.input_b{BUILD}.output_hg{OUTPUT_BUILD}.gwaslab[.added_n]
```

For example, `CAD_SCHUNKERT.EUR.input_b18.output_hg38.gwaslab`.

### Stage checkpoint files *(intermediate; safe to remove after pipeline completes)*

| File | Written by stage | Description |
|------|-----------------|-------------|
| `<stem>.preprocess.parquet` | `preprocess` | Standardised DataFrame (BROTLI-compressed) |
| `<stem>.preprocess.json` | `preprocess` | Build metadata (reference, build_num, input_build) |
| `<stem>.normalize.pkl` | `process-normalize` | After basic checks + liftover |
| `<stem>.checkref.pkl` | `process-check-ref` | After reference allele check + flip |
| `<stem>.inferstrand.pkl` | `process-infer-strand` | After 1KG strand inference |
| `<stem>.assignrsid.pkl` | `process-assign-rsid` | After dbSNP rsID assignment |

### Final output files *(never removed by cleanup)*

| File | Description |
|------|-------------|
| `<stem>.pkl` | Final raw GWASLab `Sumstats` pickle |
| `<stem>.parquet` | Raw harmonised data (parquet) |
| `<stem>.tsv.gz` | Raw harmonised data (TSV.GZ) |
| `<stem>.log` | GWASLab internal log |
| `<stem>.qc.pkl` | QC-filtered pickle |
| `<stem>.qc.parquet` | QC-filtered data (parquet) |
| `<stem>.qc.tsv.gz` | QC-filtered data (TSV.GZ) |
| `<stem>.cojo.gz` | COJO file (raw) |
| `<stem>.qc.cojo.gz` | COJO file (QC-filtered) |
| `<stem>.leads.tsv` | Lead SNPs (QC-filtered, p < 5×10⁻⁸) |
| `PLOTS/<stem>.*.png` | Diagnostic plots (Manhattan, QQ, DAF, histograms) |
| `<GWAS>.gwaslab_process.log` | Pipeline run log |

---

# 📥 Reference file management

`gwaslab.download_refs.py` is a complete inventory of all gwaslab reference files. It uses gwaslab's built-in `gl.download_ref()` function; `.tbi` index files are fetched automatically alongside VCFs. Files already present are listed as comments — only the missing 1KG population VCFs (AFR, EAS, AMR, SAS) are active by default. Uncomment any entry to (re-)download it.

## Reference file inventory

| Keyword | File | Default |
|---------|------|---------|
| `1kg_eur_hg19` / `_hg38` | `EUR.ALL.split_norm_af.*` | commented — already present |
| `1kg_pan_hg19` / `_hg38` | `PAN.ALL.split_norm_af.*` | commented — already present |
| `1kg_afr_hg19` / `_hg38` | `AFR.ALL.split_norm_af.*` | **downloaded** |
| `1kg_eas_hg19` / `_hg38` | `EAS.ALL.split_norm_af.*` | **downloaded** |
| `1kg_amr_hg19` / `_hg38` | `AMR.ALL.split_norm_af.*` | **downloaded** |
| `1kg_sas_hg19` / `_hg38` | `SAS.ALL.split_norm_af.*` | **downloaded** |
| `1kg_hm3_hg19_eaf` / `_hg38_eaf` | `PAN.hapmap3.hg{19,38}.EAF.tsv.gz` | commented — already present |
| `1kg_dbsnp151_hg19_auto` / `_hg38_auto` | `1kg_dbsnp151_hg{19,38}_auto.txt.gz` | commented — already present |
| `dbsnp_v151_hg19` / `_hg38` | `00-All.vcf.gz` (very large, NCBI FTP) | commented — already present |
| `dbsnp_v157_hg19` / `_hg38` | `GCF_000001405.{25,40}.gz` (very large, NCBI FTP) | commented — already present |
| `ucsc_genome_hg19` / `_hg38` | `hg{19,38}.fa.gz` (large, UCSC) | commented — already present |
| `recombination_hg19` / `_hg38` | `recombination_hg{19,38}.tar.gz` | commented — already present |
| `ensembl_hg19_gtf` / `_hg38_gtf` | `Homo_sapiens.GRCh3{7,8}.*.gtf.gz` | commented — already present |
| `refseq_hg19_gtf` / `_hg38_gtf` | `GRCh3{7,8}_latest_genomic.gtf.gz` | commented — already present |

> **Note:** recombination maps and GTF files are auto-downloaded by gwaslab at runtime when needed. Pre-downloading them is optional but useful on HPC nodes without outbound internet during jobs.

## Usage

```bash
conda activate gwas2cojo
python gwaslab.download_refs.py
# or with a custom path:
python gwaslab.download_refs.py --ref-dir /path/to/references/gwaslab/
```

The default target directory is `/path/to/references/gwaslab/`. Pass `--ref-dir` to override. After completion, the script prints a summary of all downloaded references via `gl.check_downloaded_ref()`.

> **Note:** Files are hosted on Dropbox and may not be accessible from all HPC networks. Run from a login node with outbound internet access, or use a node with proxy configured.

---

# 🖥️ HPC: SLURM submission

For large-scale processing of many GWAS studies on a compute cluster, four helper scripts are provided.

## ⚙️ One-time site setup

All HPC scripts read their paths and settings from a single file — **`gwas2cojo.conf`** — that lives next to the scripts. Configure it once; every script picks it up automatically.

```bash
# 1. Copy the template
cp gwas2cojo.conf.example gwas2cojo.conf

# 2. Open gwas2cojo.conf and fill in the five values:
#    PYTHON_SCRIPT  — absolute path to gwaslab.process.py
#    REF_DIR        — directory containing gwaslab reference files
#    OUT_BASE       — base output directory (per-study subdirs + SLURM logs go here)
#    CONDA_ENV      — conda environment name (default: gwas2cojo)
#    EMAIL          — your email address for SLURM failure notifications
nano gwas2cojo.conf
```

`gwas2cojo.conf` and `gwas_list.txt` are both listed in `.gitignore` so your local settings and study list are never accidentally committed. `gwas2cojo.conf.example` and `gwas_list.example.txt` (with placeholder values) are the committed templates.

## Files

| File | Description |
|------|-------------|
| `gwas2cojo.conf.example` | Site configuration template — copy to `gwas2cojo.conf` and fill in your paths |
| `gwas_list.example.txt` | Study list template (3 example studies) — copy to `gwas_list.txt` and update paths |
| `gwaslab.process.array_for_submit.sh` | SLURM worker — runs one `gwaslab.process.py` call for one study and one stage |
| `gwaslab.process.submit.sh` | Submit one full-pipeline job per study (`--stage all`) |
| `gwaslab.process.submit_staged.sh` | Submit a chained per-stage job per study with `--dependency=afterok` |
| `gwaslab.process.cleanup.sh` | Remove intermediate checkpoint files after a successful run |
| `gwaslab.download_refs.py` | Download missing 1KG population reference VCFs using `gl.download_ref()` |

## Config file format (`gwas_list.txt`)

Semicolon-separated, one GWAS per line. Lines starting with `#` and blank lines are ignored.

```
# COL1  Full path to input GWAS file (directory + filename)
# COL2  GWAS name         → --gwas and output subdirectory
# COL3  Population        → EUR / EAS / SAS / AFR / AMR / PAN
# COL4  Build             → 18 / 19 / 38
# COL5  N (total)         → integer or .
# COL6  N_cases           → integer or .
# COL7  N_controls        → integer or .
# COL8  MEM               → SLURM memory for HEAVY stages (process-infer-strand,
#                           process-assign-rsid, process-check-af); e.g. 128G or 256G
# COL9  TIME              → SLURM time limit for HEAVY stages (HH:MM:SS, e.g. 96:00:00)
# COL10 MEM_LIGHT         → SLURM memory for LIGHT stages (process-normalize,
#                           process-check-ref, qc); e.g. 32G or 64G or 128G
# COL11 TIME_LIGHT        → SLURM time limit for LIGHT stages (HH:MM:SS, e.g. 12:00:00)
```

If **all three** of COL5–COL7 are non-`.`, `--force-n` is added automatically.

> **Note on MEM_LIGHT:** Set this higher for studies with many columns or complex allele structure. Multi-ancestry meta-analyses can require 128G even at `process-check-ref` despite having fewer variants than standard EUR studies.

Example:

```
/data/gwas/CARDIoGRAM/CHD_meta_SAIGE.out.gz;CAD_Aragam;EUR;19;.;.;984168;128G;96:00:00;64G;24:00:00
/data/gwas/LIPIDS/HDL_EUR.h.tsv.gz;HDL_EUR;EUR;19;.;.;.;64G;08:00:00;32G;12:00:00
/data/gwas/AFGen/AF_TOPMed.txt.gz;AF;PAN;19;.;.;.;128G;96:00:00;128G;24:00:00
```

---

## Option 1: Single full-pipeline job per study (`gwaslab.process.submit.sh`)

Submits one SLURM job per study running `--stage all`. The `MEM` and `TIME` from the config are used for the entire job. Suitable when you want simplicity over granular resource control.

```bash
bash gwaslab.process.submit.sh gwas_list.txt
```

Extra `sbatch` arguments can be appended:

```bash
bash gwaslab.process.submit.sh gwas_list.txt --partition=highmem
```

---

## Option 2: Staged per-study chain (`gwaslab.process.submit_staged.sh`)

Submits SLURM jobs chained with `--dependency=afterok`. If a stage fails, SLURM automatically cancels all downstream stages for that study (`DependencyNeverSatisfied`). Other studies are completely independent and keep running.

### Per-chromosome pipeline (default)

The heavy VCF-sweep stages run as **SLURM array jobs** (one task per chromosome), so each job processes ~1/22 of the variants and requires proportionally less memory.

```
preprocess → process-normalize → process-split
    → [array 1-26] process-check-ref
    → [array 1-26] process-infer-strand
    → [array 1-26] process-assign-rsid  (only if --dbsnp)
    → [array 1-26] process-check-af
    → merge
```

Array task IDs: **1–22** autosomes · **23** = X · **24** = Y · **25** = nonPAR · **26** = MT.
Tasks for chromosomes absent from a study exit gracefully (exit 0), satisfying `afterok` automatically.
`afterok` on an array job ID waits for **all** 26 tasks before the next stage starts.

`merge` concatenates all per-chromosome parquets, runs QC + plots + leads + COJO, and produces the same final outputs as the old `qc` + `cojo` stages.

### Resource tiers

| Tier | Stages | Resource source |
|------|--------|----------------|
| Fixed (trivial) | `preprocess` (32G/30min), `process-split` (16G/30min) | Hardcoded in script |
| **LIGHT** | `process-normalize`, `process-check-ref` (per chr), `merge` | COL10 (`MEM_LIGHT`) + COL11 (`TIME_LIGHT`) |
| **HEAVY** | `process-infer-strand`, `process-assign-rsid`, `process-check-af` (all per chr) | COL8 (`MEM`) + COL9 (`TIME`) |

Because each per-chr job sweeps only ~1/22 of the VCF region, the HEAVY tier memory requirement is substantially lower than the equivalent whole-genome sweep.

### Submit

```bash
bash gwaslab.process.submit_staged.sh gwas_list.txt
```

Output shows the two-tier resources and the job IDs per study:

```
GWAS                    MEM_L    TIME_L      MEM_H    TIME_H      pre=... nrm=... spl=... chr=... ist=... rsi=... caf=... mrg=...
CAD_Aragam              64G      24:00:00    128G     96:00:00    pre=1001 nrm=1002 spl=1003 chr=1004 ...
```

A timestamped submission log is automatically written to `${LOG_BASE}/gwaslab.process.submit_staged_YYYYMMDD_HHMMSS.log`.

### Monitor and cancel

```bash
squeue -u $USER                                           # all jobs
scancel --name=gl_CAD_Aragam_preprocess                   # cancel one study's chain
sacct -j 1004 --format=JobID,State,Elapsed,MaxRSS        # check completed stage
```

---

## Cleanup intermediate checkpoints (`gwaslab.process.cleanup.sh`)

After a successful run the intermediate pickle checkpoints (`.normalize.pkl`, `.checkref.pkl`, `.inferstrand.pkl`, `.assignrsid.pkl`, `.preprocess.parquet`, `.preprocess.json`) can be safely removed. Final outputs (`.parquet`, `.tsv.gz`, `.qc.*`, `.cojo.gz`, `.leads.tsv`, `PLOTS/`) are **never** touched.

```bash
# Remove checkpoints for a single study
bash gwaslab.process.cleanup.sh --study CAD_Aragam

# Remove checkpoints for all studies in a config file
bash gwaslab.process.cleanup.sh --config gwas_list.txt

# Preview what would be deleted (no files are removed)
bash gwaslab.process.cleanup.sh --config gwas_list.txt --dry-run

# Also keep the final raw pickle (default: removed)
bash gwaslab.process.cleanup.sh --config gwas_list.txt --keep-raw-pkl

# Also remove the QC pickle (default: kept)
bash gwaslab.process.cleanup.sh --config gwas_list.txt --remove-qc-pkl
```

---

## 📖 License

```
The MIT License (MIT)
Copyright (c) 1979-2026 Lennart P.L. Landsmeer (lennart[at]landsmeer[dot]email), Emma J.A. Smulders (emmasmulders[at]outlook[dot]com) & Sander W. van der Laan (s.w.vanderlaan[at]gmail[dot]com).
```

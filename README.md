[gwas2cojo](https://github.com/CirculatoryHealth/gwas2cojo)<img align="right" height="200" src=logo/fulllogo_transparent.png>
============
[![Languages](https://skillicons.dev/icons?i=bash,py)](https://skillicons.dev) 

`gwas2cojo.py` is a public python script that aligns a public GWAS dataset to a genetic reference, to enable large scale cross dataset comparisons with public available GWAS datasets. Among others, it tries to deal with different dataformats and different genome builds.
Most importantly, it aligns the variant notation such that swapped, translated, wrong or ambiguous ambivalent allels are corrected or removed.

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

## ⚙️ Requirements

* Python ≥3.6
* Tested on Linux and macOS
* Optional dependencies:
* `numpy` → for allele frequency concordance (R²)
* `pyliftover` → for genome build conversion (hg18/hg19/hg38)

## 🧩 Installation

We recommend using `mamba` (a faster drop-in replacement for `conda`).

### 1️⃣ Create and activate the environment

```
mamba env create -f environment.yml
mamba activate gwas2cojo
```

### 2️⃣ Verify installation

Check that the required Python modules are available:

```
python -c "import numpy, pyliftover; print('OK')"
```

If you see OK, the environment is ready.

### 3️⃣ (Optional) Test the script

You can confirm that `gwas2cojo.py` runs correctly:

```
python gwas2cojo.py --help
```

Expected output:

A help message showing available arguments such as `--gwas`, `--gen`, `-o`, `-r`, `--fmid`, and others.

If you want to run a minimal functionality test:

```
python gwas2cojo.py --header-only --gwas example_gwas.txt.gz
```

This checks header parsing without performing alignment.

## 🔧 Troubleshooting

### 🧱 Missing `pyliftover`

If `pyliftover` is missing, install it manually:

```
mamba install -c conda-forge pyliftover
```

### 🐍 Prefer `conda` instead of `mamba`?

If you prefer `conda`:

```
conda env create -f environment.yml
conda activate gwas2cojo
```

---

# 🔬 gwaslab.process.py — GWASLab Processing Pipeline

`gwaslab.process.py` is a standalone pipeline built on the [GWASLab](https://github.com/Cloufield/gwaslab) library. It takes raw GWAS summary statistics and runs them through a fully automated processing chain: standardisation, strand inference, build liftover, dbSNP annotation, allele-frequency validation, QC filtering, and output generation in multiple formats (pickle, parquet, TSV.GZ, COJO). It is the recommended successor to `gwas2cojo.py` for new datasets.

## 🗺️ Pipeline overview

Steps run in order; individual steps can be toggled with the flags described below.

| Step | Description | Toggle |
|------|-------------|--------|
| 1 | Load and parse input file; auto-detect column names | always |
| 2 | Verify / normalise genome build | always |
| 3 | Correct P-values, SE, and N columns | always |
| 4 | Standardise to GWASLab `Sumstats` object | always |
| 5 | Plot raw input histograms | `--figures` |
| 6 | `basic_check` — flag malformed variants | always |
| 7 | `remove_dup` — drop duplicate variants | always |
| 8 | Liftover (hg18→hg38 or hg19→hg38) | `--liftover` |
| 9 | `check_ref` against hg38 FASTA | always |
| 10 | `fix_id` — normalise variant IDs | always |
| 11 | `infer_strand2` — resolve strand from 1KG VCF | always |
| 12 | `assign_rsid` — annotate with dbSNP rsIDs | `--dbsnp` |
| 13 | `check_af2` — allele-frequency concordance vs 1KG VCF | always |
| 14 | Save raw outputs (pickle + parquet + TSV.GZ) | always |
| 15 | Manhattan, QQ, and DAF plots (pre-QC) | `--figures` |
| 16 | QC filter (EAF, BETA, SE, INFO, HWE, MAC, DAF) | `--qc` |
| 17 | Save QC outputs | `--qc` |
| 18 | Manhattan, QQ, and DAF plots (post-QC) | `--figures` + `--qc` |
| 19 | Extract genome-wide significant lead SNPs | `--leads` |
| 20 | Write COJO-format file(s) | `--cojo` |

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

### Resume from a saved pickle (skip re-processing)

If the pipeline completed successfully up to the pickle save step but failed later (e.g. during plotting), you can reload from the pickle and re-run downstream steps without repeating the expensive processing:

```bash
python3 gwaslab.process.py \
    --gwas    CAD_SCHUNKERT \
    --input   cardiogram_gwas_results_edited.txt.gz \   # ignored; pickle is used instead
    --directory /data/gwas/CARDIoGRAM \                # ignored; pickle is used instead
    --ref     /data/references/gwaslab \
    --output  /data/results/CAD_SCHUNKERT \            # must match the original run
    --population EUR \
    --build   18 \                                     # must match the original run
    --liftover \                                       # must match the original run
    --qc --figures --leads \
    --only-qc
```

> **Important:** `--output`, `--population`, `--build`, and `--liftover` must match the original run exactly — they determine the pickle filename that is looked up.

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
| `--only-qc` | Skip loading; reload from an existing pickle and run QC / plots / leads only |
| `--fill-eaf` | Look up missing EAF values from the 1KG reference VCF |
| `--figures` | Generate diagnostic plots (Manhattan, QQ, DAF, histograms) |
| `--leads` | Extract genome-wide significant lead SNPs (p < 5×10⁻⁸) |

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

| File | Description |
|------|-------------|
| `<stem>.pkl` | GWASLab `Sumstats` pickle — used for `--only-qc` resume |
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

# 🖥️ HPC: SLURM array-job submission

For large-scale processing of many GWAS studies on a compute cluster, two helper scripts are provided.

## Files

| File | Description |
|------|-------------|
| `gwaslab_array.sh` | SLURM batch script; runs one `gwaslab.process.py` call per array task |
| `submit_gwaslab.sh` | Submission helper; counts config entries and calls `sbatch --array=1-N` |
| `gwas_list.txt.example` | Annotated example config file |

## Config file format

Tab-separated, one GWAS per line. Lines starting with `#` and blank lines are ignored.

```
# COL1  Full path to input GWAS file
# COL2  GWAS name          → --gwas and output subdirectory
# COL3  Population         → EUR / EAS / SAS / AFR / AMR / META
# COL4  Build              → 18 / 19 / 38
# COL5  N total            → integer, or . if not applicable
# COL6  N cases            → integer, or . if not applicable
# COL7  N controls         → integer, or . if not applicable
```

If **all three** of COL5–COL7 are non-`.`, `--force-n` is added automatically.

Example:

```
/data/gwas/CARDIoGRAM/cardiogram_gwas.txt.gz	CAD_SCHUNKERT	EUR	18	.	.	.
/data/gwas/MI/nikpay2015.txt.gz	MI_NIKPAY	EUR	19	.	43676	128199
/data/gwas/T2D/mahajan2018.txt.gz	T2D_MAHAJAN	EUR	38	898130	.	.
```

## Quick start

**1.** Copy and edit the example config:

```bash
cp gwas_list.txt.example gwas_list.txt
# fill in your actual paths and study details
```

**2.** Edit the `USER CONFIGURATION` block at the top of `gwaslab_array.sh` for your cluster:

```bash
PYTHON_SCRIPT="/hpc/...path.../gwaslab.process.py"
REF_DIR="/hpc/...path.../references/gwaslab"
OUT_BASE="/hpc/...path.../results"
CONDA_ENV="gwas2cojo"
```

**3.** Submit:

```bash
bash submit_gwaslab.sh gwas_list.txt
```

This counts valid entries in the config and submits one SLURM array task per GWAS study. Extra `sbatch` arguments can be passed to override defaults:

```bash
bash submit_gwaslab.sh gwas_list.txt --partition=highmem --time=48:00:00
```

## SLURM defaults

The defaults set in `gwaslab_array.sh` are:

| Resource | Default |
|----------|---------|
| Memory | 64 GB |
| CPUs | 4 |
| Wall time | 24 h |

Adjust the `#SBATCH` directives in `gwaslab_array.sh` or override on the command line when submitting.

---

## 📖 License

```
The MIT License (MIT)
Copyright (c) 1979-2026 Lennart P.L. Landsmeer (lennart[at]landsmeer[dot]email), Emma J.A. Smulders (emmasmulders[at]outlook[dot]com) & Sander W. van der Laan (s.w.vanderlaan[at]gmail[dot]com).
```

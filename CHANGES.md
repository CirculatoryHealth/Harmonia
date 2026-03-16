# CHANGES

This document tracks changes to the codebase. Each entry should include a brief description of the change, the files affected, and any relevant context or reasoning behind the change. This helps maintain a clear history of modifications and facilitates collaboration among developers.

## 2026-03-16 🆕 Fine-grained process sub-stages, staged submit, and cleanup (v1.3.0)
- 🛠️**Updated**: Bumped version to `1.3.0` (`2026-03-16`).
- 🆕**Added**: Five `--stage process-*` sub-stages to `gwaslab.process.py`, splitting the monolithic process stage by memory profile. Each sub-stage saves a pickle checkpoint so subsequent stages can be submitted as independent SLURM jobs with their own resources:
    - `process-normalize`    — `basic_check` + `remove_dup` + `liftover` → `{stem}.normalize.pkl` (medium)
    - `process-check-ref`    — `check_ref` + `flip_allele_stats` + `fix_id` → `{stem}.checkref.pkl` (medium)
    - `process-infer-strand` — `infer_strand2` + `flip_allele_stats` → `{stem}.inferstrand.pkl` (high — 1KG VCF sweep)
    - `process-assign-rsid`  — `assign_rsid` via dbSNP VCF sweep → `{stem}.assignrsid.pkl` (extreme — dbSNP sweep; skipped when `--dbsnp` not set)
    - `process-check-af`     — `check_af2` → final raw outputs `.pkl` + `.parquet` + `.tsv.gz` (high — 1KG VCF sweep)
- 🆕**Added**: `run_normalize()`, `run_check_ref()`, `run_infer_strand()`, `run_assign_rsid()`, `run_check_af()` — individual runner functions extracted from `run_processing()`, each containing exactly the steps for their sub-stage.
- 🆕**Added**: `save_process_checkpoint()` and `load_process_checkpoint()` — pickle-based checkpoint I/O for process sub-stages, with descriptive error messages on missing files.
- 🆕**Added**: `_PROCESS_CHECKPOINT_META` lookup dict mapping checkpoint suffixes to stage names and their predecessor, used in error messages when a checkpoint is missing.
- 🛠️**Updated**: `--stage` choices in `parse_args()` now include all five `process-*` sub-stages alongside `all`, `preprocess`, `qc`, and `cojo`.
- 🛠️**Updated**: `main()` — added guard that exits with an error if `--stage process-assign-rsid` is used without `--dbsnp`.
- 🛠️**Updated**: `main()` — `_pickle_required_stages` guard extended to cover all process sub-stages that require a prior-stage checkpoint.
- 🛠️**Updated**: Header comment in `main()` documents the full checkpoint chain: `preprocess → normalize → checkref → inferstrand → assignrsid → checkaf → qc → cojo`.
- 🆕**Added**: `gwaslab.process.submit_staged.sh` — new script that submits one SLURM job per stage per study with `--dependency=afterok` chaining. If a stage fails, SLURM cancels all downstream stages for that study automatically; other studies are unaffected. `MEM` and `TIME` from `gwas_list.txt` are applied to the two heaviest stages (`process-infer-strand` and `process-assign-rsid`); all other stages use fixed resource defaults defined at the top of the script. `process-assign-rsid` is omitted when `--dbsnp` is absent from `WORKER_FLAGS`.
- 🆕**Added**: `gwaslab.process.cleanup.sh` — removes all intermediate checkpoint files after a successful run. Final outputs (`.parquet`, `.tsv.gz`, `.qc.*`, `.cojo.gz`, `.leads.tsv`, `.log`, `PLOTS/`) are never touched. Supports `--study NAME`, `--all`, or `--config gwas_list.txt` scope; `--dry-run` prints what would be deleted without removing; `--keep-raw-pkl` preserves the final raw pickle; `--remove-qc-pkl` also removes the QC pickle (kept by default).

## 2026-03-16 🆕 Pipeline staging in gwaslab.process.py (v1.2.0)
- 🛠️**Updated**: Bumped version to `1.2.0` (`2026-03-16`).
- 🆕**Added**: `--stage` flag to `gwaslab.process.py` with four stages: `preprocess`, `process`, `qc`, and `cojo` (plus `all`, the default, which preserves the existing end-to-end behaviour). Each stage can be submitted as a separate SLURM job with its own `--mem` and `--time`, allowing resource-light stages to run at `64G / 48h` while memory-intensive steps (`process`: `check_ref`, `infer_strand2`, `assign_rsid`) can be given `128G–256G / 96h` independently.
- 🆕**Added**: `save_preprocess_checkpoint()` — writes `{stem}.preprocess.parquet` (BROTLI-compressed standardised DataFrame) and `{stem}.preprocess.json` (detected build metadata) as the handoff from `--stage preprocess` to `--stage process`.
- 🆕**Added**: `load_preprocess_checkpoint()` — reads the parquet + JSON checkpoint written by `--stage preprocess` and restores the `reference`, `build_num`, and `input_build` so subsequent stages use identical file stems.
- 🆕**Added**: `import json` to top-level imports (previously absent; required by the new checkpoint metadata functions).
- 🛠️**Updated**: `main()` — refactored into four clearly labelled stage blocks (`STAGE: preprocess`, `STAGE: process`, `STAGE: qc`, `STAGE: cojo`). In `--stage all` mode the blocks execute in sequence without touching disk checkpoints, preserving current behaviour. In individual-stage mode each block saves its checkpoint and returns early.
- 🛠️**Updated**: `main()` — `--only-qc` is now a backward-compatible alias for `--stage qc`; a deprecation notice is logged when it is used.
- 🆕**Added**: Guard in `main()` that exits with an error if `--stage qc` or `--stage cojo` is combined with `--no-pickle`, since both stages require a pickle written by a prior stage.
- 🛠️**Updated**: Stage summary at end of each stage block now logs the next recommended stage invocation (e.g. `Next: --stage process  (pass the same --gwas / --build / --liftover / --output flags)`).

## 2026-03-16 🛠️ Memory efficiency improvements in gwaslab.process.py (v1.1.0)
- 🛠️**Updated**: Bumped version to `1.1.0` (`2026-03-16`).
- 🆕**Added**: `import gc` (previously commented out) to enable explicit garbage collection at stage boundaries.
- 🆕**Added**: `--no-pickle` flag to `gwaslab.process.py`. When set, `.pkl` files are skipped for both raw and QC outputs, reducing peak memory and disk usage on the save step. The gwaslab `.log` file is still written regardless. Note: `--only-qc` requires a pickle from a prior run, so it is incompatible with `--no-pickle`.
- 🛠️**Updated**: `save_raw_outputs()` and `save_qc_outputs()` — eliminated the parquet read-back pattern (`pd.read_parquet(parquet_path)`) that was used to generate the TSV.GZ. Both functions now write the TSV directly from the in-memory `gwas_obj.data` / `gwas_obj_qc.data`, avoiding a full extra copy of the data just to write one file.
- 🛠️**Updated**: `main()` — added `del gwas_data; gc.collect()` immediately after `plot_raw_histograms()` (the last use of the raw DataFrame). This frees the raw pandas DataFrame before the heavy processing steps (`check_ref`, `infer_strand2`, `assign_rsid`, `check_af2`), preventing two full-size DataFrames from coexisting in RAM throughout the pipeline.
- 🛠️**Updated**: `main()` — added `del gwas_obj; gc.collect()` immediately after `apply_qc()` returns `gwas_obj_qc`. The unfiltered object is freed before saving QC outputs and generating QC plots, so only one copy of the data is in memory at a time during the QC stage.
- 🛠️**Updated**: `write_cojo()` — removed unnecessary `.copy()` call (`df = gwas_obj.data.copy()` → `df = gwas_obj.data`). All downstream accesses are read-only (column selection, `astype`, constructing a new `pd.DataFrame`), so the copy was wasted memory.

## 2026-03-15 🛠️ Overhaul of SLURM submission and GWAS list
- 🛠️**Updated**: The `gwas_list.txt` file to use semicolons (`;`) as the field delimiter instead of tabs, avoiding parsing issues when paths or values contain whitespace.
- 🆕**Added**: Two new columns to `gwas_list.txt`: `MEM` (COL8, SLURM memory per job, e.g. `64G` or `128G`) and `TIME` (COL9, SLURM time limit per job, e.g. `48:00:00`), allowing resource requirements to be set individually per dataset.
- 🛠️**Updated**: `gwaslab.process.submit.sh` to submit one independent SLURM job per dataset instead of a single array job. Memory (`--mem`) and time (`--time`) are now read from the config file and passed to each `sbatch` call individually, so datasets with different resource needs no longer share a single limit. Each job receives its own `--job-name`, `--output`, and `--error` derived from the dataset name.
- 🛠️**Updated**: `gwaslab.process.array_for_submit.sh` to act as a single-dataset worker script. Removed array job logic (`SLURM_ARRAY_TASK_ID`), removed fixed `--mem`, `--time`, `--output`, and `--error` SBATCH directives (these are now set dynamically by `gwaslab.process.submit.sh`). The script now accepts a semicolon-delimited config line as its first argument and parses it directly.

## 2025-03-12 🛠️ Updates to GWAS list
- 🆕**Added**: New GWAS datasets to the `gwas_list.txt` file, including:
    - AFGen Roselli 2018 dataset for allele frequencies (AF) with b38 positions.
    - GLGC Graham 2021 datasets for HDL, LDL, TC, TG, and non-HDL traits in European populations.
- 🧰**Fixed**: Issue with time of the SLURM job in `gwaslab.process.array_for_submit.sh` to allow for longer processing times, especially for larger GWAS datasets. Updated the time limit from 1 hour to 4 hours to accommodate the increased computational demands of processing multiple large GWAS datasets.

## 2025-03-12 🛠️ Updates to GWAS list
- 🆕**Added**: New GWAS datasets to the `gwas_list.txt` file, including:
    - ISGC GigaStroke datasets for ALLSTROKE, IS, CES, LAS, and SVD subtypes.
    - CHARGE cIMT (Franceschini 2018) and CHARGE Plaque (Franceschini 2018) datasets.
- 🛠️**Updated**: The `gwas_list.txt` file to ensure consistency in formatting and correct file paths.
- 🛠️**Updated**: Changed the SLURM parameters for `gwaslab.process.array_for_submit.sh`. 

## 2025-03-12 🆕 New functions
- 🆕**Added**: A notebook to test drive some functions and option using `gwaslab`. 
    - New functionality to save QC-filtered output in `gwaslab.process.ipynb`.
    - Plots for QC-filtered dataset in `gwaslab.process.ipynb`.
    - Extraction of lead SNPs in `gwaslab.process.ipynb`.
- 🆕**Added**: New script to process a given GWAS using `gwaslab`. This script will:
    - Load the GWAS summary statistics.
    - Perform liftover if necessary.
    - Check reference alleles and flip if needed.
    - Check for duplicates and remove them.
    - Check for strand issues and resolve them.
    - Check for allele frequency issues and filter variants accordingly.
    - Perform QC filtering.
    - Generate plots for both the full dataset and the QC-filtered dataset.
    - Extract lead SNPs from the QC-filtered dataset.
    - Ensure the `stem` variable is defined for both normal and --only-qc paths, allowing consistent file naming across different branches of the code.
    - Updated plotting functions in `gwaslab.process.py` to include verbose logging and ensure that plots are saved with the correct DPI settings.
    - Handles the case where a pickle file was created and the --only-qc flag is used to regenerate plots without re-running the full pipeline.
- 🛠️**Updated**: The `LICENSE` file to correct the copyright year.
- 🛠️**Updated**: The `.gitignore` file to include new directories and files that should be ignored by git.
- 🛠️**Updated**: The `CHANGES.md` file to document the new functions and updates made to the codebase.
- 🛠️**Updated**: The `README.md` file to reflect the new functionality and provide instructions for using the new script and notebook.
- 🛠️**Updated**: The `gwaslab.process.py` file to include the new script for processing GWAS summary statistics and to ensure that the `stem` variable is defined in all relevant branches of the code.
- 🆕**Added**: Scripts for submitting GWAS processing jobs:
    - `gwaslab.process.submit.sh`: A shell script to submit a GWAS processing job to a cluster using `sbatch`.
    - `gwaslab.process.array_for_submit.sh`: A shell script to submit an array of GWAS processing jobs for multiple datasets or parameters. This is controlled by the `gwaslab.process.submit.sh` script, which can be configured to run multiple instances of the processing script with different arguments.
    - `gwas_list.txt`: A text file containing a list of GWAS datasets to be processed. This file is used by the `gwaslab.process.array_for_submit.sh` script to determine which datasets to process in the array job. Each line in the file should specify a GWAS dataset, and the processing script will read this file to know which datasets to run on.

# CHANGES

This document tracks changes to the codebase. Each entry should include a brief description of the change, the files affected, and any relevant context or reasoning behind the change. This helps maintain a clear history of modifications and facilitates collaboration among developers.

## 2026-03-18 🐛 SIGILL crash on older HPC compute nodes — polars requires AVX2 (v1.4.6)
- **Bug**: All per-chromosome SLURM array jobs crashed with `Illegal instruction (core dumped)` preceded by a `polars` RuntimeWarning: `Missing required CPU features: avx2, fma, bmi1, bmi2, lzcnt, movbe`. The standard `polars` wheel on PyPI is compiled with AVX2/FMA intrinsics; older HPC compute nodes (pre-Haswell microarchitecture) lack those instruction-set extensions. Setting `POLARS_SKIP_CPU_CHECK=1` only suppresses the Python-level RuntimeWarning — the binary still faults the moment any AVX2 instruction executes. The crash occurs inside `gwaslab` which imports `polars` internally.
- 🐛**Fixed** `environment.yml` — replaced `"polars>=1.27.0"` with `"polars-lts-cpu>=1.27.0"`. The `polars-lts-cpu` PyPI package is the official CPU-compatible build of Polars, compiled for SSE2/SSE4 without AVX2 or FMA requirements. It provides an identical public API and can run on any x86-64 node regardless of CPU generation. The tradeoff is a modest performance reduction on modern nodes (typically 10–20 % slower for vectorised operations), which is negligible compared to the I/O and Python overhead in this pipeline.
- ⚠️ **Action required**: rebuild the conda environment after pulling this change: `mamba env remove -n gwas2cojo && mamba env create -f environment.yml`. If the environment must be updated in-place without rebuilding: `pip uninstall polars && pip install "polars-lts-cpu>=1.27.0"`.
- 🛠️**Updated**: Bumped version to `1.4.6` (`2026-03-18`).

## 2026-03-18 🐛 Categorical re-encoding by gwaslab __init__ breaks flip_allele_stats in per-chr stages (v1.4.5)
- **Bug**: `process-check-ref` (and downstream per-chromosome stages) crashed with `TypeError: Cannot setitem on a Categorical with a new category, set the categories first` inside `gwaslab.flip_allele_stats()`, even though `load_chrom_parquet()` already converted Categorical columns to `object` before passing the DataFrame to `make_sumstats_from_chrom_df()`. Root cause: gwaslab's `gl.Sumstats.__init__()` calls `basic_check()` internally, which re-encodes `EA` and `NEA` as `pd.Categorical`. Because the input is a per-chromosome shard, each column's category set only contains the allele values observed on that chromosome. When `flip_allele_stats()` then tries to swap `EA↔NEA` for 117,715 variants (e.g. LAS chromosome 9), it attempts to assign an `EA` value that is present in `NEA`'s category set but absent from `EA`'s subset, causing the pandas error. The earlier `v1.4.2` fix in `load_chrom_parquet()` was not sufficient because it converted before construction, and construction undoes it.
- 🐛**Fixed** `gwaslab.process.py` — `make_sumstats_from_chrom_df()` now runs a second Categorical-to-object conversion on `gwas_obj.data` immediately *after* `make_sumstats_object()` returns. All columns identified by `select_dtypes(include="category")` (typically `EA`, `NEA`, `SNPID`) are converted to plain `object` dtype. This conversion is applied once at object-creation time and persists through all downstream per-chromosome processing steps (`check_ref`, `infer_strand`, `assign_rsid`, `check_af`). The `load_chrom_parquet()` conversion is kept as a pre-construction safety net.
- 🛠️**Updated**: Bumped version to `1.4.5` (`2026-03-18`).

## 2026-03-18 🧪 Synthetic test datasets for gwas2cojo.py
- 🆕**Added** `test/generate_test_data.py` — stdlib-only Python script that generates 10,000 synthetic biallelic SNPs and writes four files: a genetic reference (`ref.txt.gz`) and three GWAS summary-statistics files covering different column-naming conventions (`gwas_metal_tab.txt.gz`, `gwas_plink2.txt.gz`, `gwas_saige.txt.gz`). Datasets include deliberate complement-strand and allele-switched variants so the NOP, FLIP, and translated-allele branches of `select_action()` are all exercised. Regenerate at any time with `python3 test/generate_test_data.py`.
- 🆕**Added** `test/ref.txt.gz`, `test/gwas_metal_tab.txt.gz`, `test/gwas_plink2.txt.gz`, `test/gwas_saige.txt.gz` — pre-generated test fixtures (~870 KB total) committed so `bash test/run.sh` runs without a generation step.
- 🛠️**Updated** `test/run.sh` — extended from a 1-variant smoke test to four tests: (0) the original A.txt/B.txt sanity check; (1) METAL-style with full auto-detection and `gwas2cojo-verify.py` validation (0 errors); (2) PLINK2-style exercising the new `A1FREQ`/`OBS_CT` aliases; (3) SAIGE-style requiring `--gwas:effect Allele2 --gwas:other Allele1 --gwas:freq AF_Allele2` manual overrides.

## 2026-03-18 🐛 Bug fixes and column alias alignment in gwas2cojo.py and gwas2cojo-verify.py

### `gwas2cojo-verify.py` → `v1.0.1` — logic errors in allele-comparison helpers
- 🐛**Fixed** `equal_alleles(a, b)` — the second comparison was `a.ref == b.ref` (duplicate of the first), meaning the *other* allele (`oth`) was never checked. Corrected to `a.oth == b.oth`.
- 🐛**Fixed** `switched_alleles(a, b)` — the body referenced free variables `gen` and `gwas` (local names inside the caller `verify()`) instead of the function parameters `a` and `b`. As a module-level function, Python resolves free variables in the global scope, so any call that reached a FLIP assertion would raise `NameError: name 'gen' is not defined`. Corrected to `a.ref == b.oth and a.oth == b.ref`.
- 🛠️**Updated** `Last update` date to `2026-03-18`.

### `gwas2cojo.py` → `v1.4.4` — four bug fixes and column alias improvements
- 🐛**Fixed** `select()` inside `read_gwas()` — `except IndexError` was catching the wrong exception type: `list.index()` raises `ValueError`, not `IndexError`. A bad user-supplied `--gwas:<col>` value therefore propagated as an uncaught `ValueError` instead of printing the helpful "not found" diagnostic. Corrected to `except ValueError`.
- 🐛**Fixed** `GWAS_H_NCONTROL_OPTIONS` / `GWAS_H_NCASE_OPTIONS` — the entries `'TotalCases'` and `'TotalSampleSize'` were in the *wrong lists* (swapped). `TotalCases` belongs in the case count list; `TotalSampleSize` belongs in the total-N list. Corrected: `GWAS_H_NCONTROL_OPTIONS` now contains `'TotalControls'`/`'n_controls'`; `GWAS_H_NCASE_OPTIONS` now contains `'TotalCases'`/`'n_cases'`.
- 🐛**Fixed** `gwas_header_auto(gwas_filename)` — function body used the undefined name `filename` instead of the parameter `gwas_filename`, and called `fopen(filename, 'rt')` with two arguments while `fopen()` only accepts one. Also used the undefined name `headers` instead of `header`. Corrected to `fopen(gwas_filename)` and `len(header)`.
- 🆕**Added** column aliases for widely-used GWAS tool outputs:
  - `GWAS_H_FREQ_OPTIONS`: `'A1FREQ'` (PLINK2 `.afreq`/`.linear`/`.logistic`), `'FRQ'` (PLINK 1.9 `.frq`)
  - `GWAS_H_NTOTAL_OPTIONS`: `'OBS_CT'` (PLINK2 observation count), `'n_total'`
- 🛠️**Updated** `Last update` date to `2026-03-18`.

## 2026-03-18 🐛 gwas2cojo.conf not found when running as SLURM job
- 🛠️**Updated**: Bumped version to `1.4.4` (`2026-03-18`).
- **Bug**: After the `gwas2cojo.conf` introduction, SLURM jobs immediately failed with `ERROR: /var/spool/slurmd/job<ID>/gwas2cojo.conf not found`. SLURM copies the worker script (`array_for_submit.sh`) to its own temporary spool directory before executing it on the compute node, so `BASH_SOURCE[0]` inside the job resolves to the spool path rather than the original script location. The conf-file lookup `"${SCRIPT_DIR}/gwas2cojo.conf"` therefore searched in `/var/spool/slurmd/job<ID>/` where no conf file exists.
- 🐛**Fixed** `gwaslab.process.array_for_submit.sh` — the conf-loading stanza now checks the environment variable `GWAS2COJO_CONF` first (exported by the submit scripts, which run on the login node and always have the correct absolute path). The `BASH_SOURCE`-relative lookup is retained as a fallback for direct local invocation only. An improved error message names all three possible causes when the conf is still not found.
- 🛠️**Updated** `gwaslab.process.submit_staged.sh`, `gwaslab.process.submit.sh` — both scripts now `export GWAS2COJO_CONF="${CONF}"` immediately after sourcing the conf. SLURM propagates all exported environment variables to job environments by default (`--export=ALL`), so the absolute path is reliably available inside every job regardless of which spool directory SLURM uses.

## 2026-03-18 🔒 Removed hardcoded site-specific paths; added gwas2cojo.conf
- 🛠️**Updated**: Bumped version to `1.4.3` (`2026-03-18`).
- 🔒**Removed** all hardcoded HPC-specific paths and institutional email addresses (`@umcutrecht.nl`) from every tracked file in the repository so the codebase is clean for public use.
- 🆕**Added** `gwas2cojo.conf.example` — a single site-configuration template containing five variables (`PYTHON_SCRIPT`, `REF_DIR`, `OUT_BASE`, `CONDA_ENV`, `EMAIL`). Users copy it to `gwas2cojo.conf` (gitignored) and fill in their local values once.
- 🛠️**Updated** `gwaslab.process.array_for_submit.sh`, `gwaslab.process.submit.sh`, `gwaslab.process.submit_staged.sh`, `gwaslab.process.cleanup.sh` — replaced per-script `USER CONFIGURATION` blocks with a uniform conf-loading stanza (`source "${SCRIPT_DIR}/gwas2cojo.conf"`). All four scripts now emit a clear error with instructions if `gwas2cojo.conf` is missing.
- 🛠️**Updated** `.gitignore` — added `gwas2cojo.conf` and `gwas_list.txt` so local site settings and study lists are never accidentally committed.
- 🛠️**Updated** `gwaslab.download_refs.py` — replaced hardcoded `DEFAULT_REF_DIR` path and docstring with placeholder values.
- 🆕**Added** `gwas_list.example.txt` — a minimal three-study example config (`CAD_Aragam`, `CHARGE_CAC_EA`, `AF`) with placeholder `/path/to/gwas_datasets/` prefixes, HEADER comments, and resource annotations. Serves as the committed template; users copy to `gwas_list.txt` (gitignored) and update paths.
- 🛠️**Updated** `gwas2cojo.py`, `gwas2cojo-verify.py` — replaced `@umcutrecht.nl` banner addresses with obfuscated personal addresses (`lennart[at]landsmeer[dot]email`, `s.w.vanderlaan[at]gmail[dot]com`), matching the format already used in `gwaslab.process.py` and the README licence block.
- 🛠️**Updated** `README.md` — added a `⚙️ One-time site setup` section explaining the `gwas2cojo.conf` workflow; added `gwas2cojo.conf.example` and `gwas_list.example.txt` to the HPC files table; replaced remaining HPC paths and institutional emails throughout.

## 2026-03-18 🐛 Categorical dtype crash in per-chr stages and ZeroDivisionError in plot_daf
- 🛠️**Updated**: Bumped version to `1.4.2` (`2026-03-18`).
- **Bug 1**: All per-chromosome `process-check-ref` and `process-infer-strand` jobs failed with `TypeError: Cannot setitem on a Categorical with a new category, set the categories first` inside gwaslab's `flip_allele_stats()`. gwaslab encodes EA/NEA/SNPID as `pd.Categorical` after `basic_check()` for memory efficiency, and parquet round-trips preserve that dtype. A per-chromosome shard's Categorical column only contains the allele categories actually present on that chromosome. When `flip_allele_stats` tries to swap EA↔NEA for a variant whose allele (e.g. an indel sequence) exists in EA's category set but not NEA's on that chromosome, pandas refuses the assignment. In the whole-genome path this never surfaced because the full Categorical across all chromosomes includes all values in both columns simultaneously.
- 🐛**Fixed** `gwaslab.process.py` — `load_chrom_parquet()` now converts all `pd.CategoricalDtype` columns to plain `object` dtype immediately after reading the parquet (`df.select_dtypes(include="category")`), before the DataFrame is passed to `make_sumstats_from_chrom_df()`. This only affects EA/NEA/SNPID-style string columns; numeric columns (STATUS `Int64`, CHR, POS, BETA, SE, P, N, etc.) are not Categorical and are completely unaffected. gwaslab operates identically on `object` dtype allele strings for all harmonise/check operations; Categorical is purely a memory optimisation that is not required for correctness.
- **Bug 2**: The `merge` stage failed for CHARGE_CAC_EA_AA with `ZeroDivisionError: division by zero` inside gwaslab's `plot_daf()` (`num / len(sumstats)` where `len(sumstats) == 0`). The study had too few variants with a valid DAF value after processing (EAF largely absent or all-NaN), leaving an empty subset after DAF filtering inside gwaslab's plot routine.
- 🐛**Fixed** `gwaslab.process.py` — wrapped `gwas_obj.plot_daf()` in both `plot_full_dataset()` and `plot_qc_dataset()` with `try/except ZeroDivisionError`. When triggered, a `logging.warning` is emitted and the DAF plot is skipped; the rest of the merge stage (Manhattan, QQ, QC, leads, COJO) continues normally.
- 🛠️**Updated** `gwaslab.process.submit.sh` — added `NODES`, `CPUS`, `EMAIL`, and `MAIL_TYPE` variables and passed `--nodes`, `--cpus-per-task`, `--mail-type`, `--mail-user` to the `sbatch` call. Previously these were absent, relying on the now-removed `#SBATCH` directives in `array_for_submit.sh`.
- 🛠️**Updated** `gwaslab.process.array_for_submit.sh` — removed `#SBATCH --mail-type=END,FAIL` and `#SBATCH --mail-user` directives from the worker script header. SLURM merges `#SBATCH` directives from the script with command-line flags rather than letting the command line override them, so the hardcoded `END` in the script was causing end-of-job emails despite both submit scripts setting `--mail-type=FAIL`. Mail settings are now solely controlled by the calling submit script. Updated the comment to correctly name both `submit.sh` and `submit_staged.sh` as the controlling scripts.

## 2026-03-18 🐛 Missing --nodes, --cpus-per-task, and --mail-* in gwaslab.process.submit_staged.sh
- 🛠️**Updated**: Bumped version to `1.4.1` (`2026-03-18`).
- **Bug**: The v1.4.0 per-chromosome refactor of `gwaslab.process.submit_staged.sh` dropped three SLURM job settings that were present in the earlier script: `--nodes`, `--cpus-per-task`, and `--mail-type`/`--mail-user`. As a result all submitted jobs would inherit SLURM defaults (typically 1 CPU, which starves the multi-threaded Python worker that requests `--threads 8` via `WORKER_FLAGS`), and no failure-notification emails would be sent.
- 🐛**Fixed** `gwaslab.process.submit_staged.sh` — added four variables to the USER CONFIGURATION block (`NODES=1`, `CPUS=8`, `EMAIL`, `MAIL_TYPE="FAIL"`) and passed `--nodes`, `--cpus-per-task`, `--mail-type`, `--mail-user` to all eight `sbatch` calls (preprocess, normalize, split, check-ref, infer-strand, assign-rsid, check-af, merge). `CPUS` is intentionally kept in sync with the `--threads N` value in `WORKER_FLAGS`.

## 2026-03-18 🐛 Wrong build passed to reference-checking stages after liftover
- **Bug**: In all staged pipeline paths (`process-check-ref`, `process-infer-strand`, `process-assign-rsid`, `process-check-af`, `qc`, and the new `merge`), `normalise_build(REFERENCE)` was used instead of `build_num` when selecting the reference FASTA, dbSNP VCF, and Sumstats build, and when setting the chromosome map for Manhattan/QQ plots. `REFERENCE` is set from `args.build` (the original input build, e.g. `"19"`) and is never updated between staged invocations. `build_num` is correctly set to `"38"` at startup for any hg19/hg18+liftover study.
- **Impact**: For any study submitted with `--liftover`, all four heavy stages and both plot-generating stages would use:
  - `hg19.fa.gz` instead of `hg38.fa.gz` in `check_ref` → coordinates are hg38, FASTA is hg19 → nearly all variants incorrectly flagged as MISREF and lost.
  - `GCF_000001405.25.gz` (hg19 dbSNP) instead of `GCF_000001405.40.gz` (hg38 dbSNP) in `assign_rsid` → rsIDs assigned from the wrong coordinate space.
  - `build="19"` in reconstructed `Sumstats` objects (per-chr and merge paths) → wrong internal build attribute for all downstream gwaslab operations.
  - `build=reference` in `plot_mqq` → hg19 chromosome-length map applied to hg38 positions → distorted Manhattan plots.
- **Why not seen before**: The staged whole-genome path (pre-v1.4.0) always OOM'd inside `process-check-ref` or later, so these stages never produced output. The per-chromosome refactor (v1.4.0) is specifically designed to make these stages complete — meaning the wrong results would be written and stored for the first time.
- 🐛**Fixed** `gwaslab.process.py` — replaced `normalise_build(REFERENCE)` / `REFERENCE` with `build_num` at nine call sites across `main()` and `run_merge()`:
  - `process-check-ref` per-chr and whole-genome: `run_check_ref(gwas_obj, build_num, args.ref)` (FASTA path)
  - `process-assign-rsid` per-chr and whole-genome: `run_assign_rsid(gwas_obj, build_num, args.ref, …)` (dbSNP VCF path)
  - `make_sumstats_from_chrom_df(df, build_num)` in all four per-chr stage branches and in `run_merge()`
  - `plot_full_dataset(…, build_num, …)` and `plot_qc_dataset(…, build_num, …)` in `process-check-af` (whole-genome), `qc`, and `merge`

## 2026-03-17 🆕 Per-chromosome array-job pipeline for heavy stages (v1.4.0)
- 🛠️**Updated**: Bumped version to `1.4.0` (`2026-03-17`).
- **Context**: Studies were OOM-failing at `process-check-ref`, `process-infer-strand`, `process-assign-rsid`, and `process-check-af` even at 128–256 G. The root cause is that these stages sweep large VCF files (1KG ~84 M variants, dbSNP ~1 B variants) against the full genome-wide GWAS dataset. The fix splits the dataset by chromosome before the heavy stages so each VCF-sweep job works on ~1/22 of the variants.
- 🆕**Added** `process-split` stage to `gwaslab.process.py` — loads `{stem}.normalize.pkl`, splits by chromosome into per-chromosome BROTLI parquets (`{stem}.chr{N}.normalize.parquet`, N = 1–26), and writes a `{stem}.chrsplit.json` manifest. CHR values follow gwaslab's `Int64` convention: 1–22 = autosomes, 23 = X, 24 = Y, 25 = nonPAR, 26 = MT. Parquets preserve the STATUS bitmask column so gwaslab state is maintained across the per-chr jobs.
- 🆕**Added** `merge` stage to `gwaslab.process.py` — concatenates all `{stem}.chr{N}.checkaf.parquet` shards into a single genome-wide DataFrame, recreates a gwaslab `Sumstats` object (with STATUS restored), then runs QC filtering, plots (Manhattan, QQ, DAF), lead-variant extraction, and COJO output. Replaces the separate `qc` + `cojo` stages in the per-chromosome pipeline path.
- 🆕**Added** `--chrom N` argument (int 1–26) to `gwaslab.process.py`. When set, `process-check-ref`, `process-infer-strand`, `process-assign-rsid`, and `process-check-af` each operate on a single chromosome shard (`{stem}.chr{N}.{prev}.parquet` → `{stem}.chr{N}.{next}.parquet`). If the shard does not exist the stage exits gracefully with exit code 0, satisfying SLURM `afterok` dependencies for the next array stage.
- 🆕**Added** helper functions: `load_chrom_parquet()`, `save_chrom_parquet()`, `make_sumstats_from_chrom_df()`, `split_by_chrom()`, `load_chrsplit_manifest()`, `run_merge()`.
- 🛠️**Updated** `gwaslab.process.submit_staged.sh` — the four heavy process stages are now submitted as SLURM array jobs (`--array=1-26`); a `process-split` job is inserted between `process-normalize` and the array stages; a `merge` job replaces the `qc` + `cojo` tail. Fixed resources: `process-split` 16 G / 30 min. `afterok` on an array job ID waits for all 26 tasks; absent-chromosome tasks (exit 0) satisfy the dependency automatically. Job count per study: ~107 (vs. 8 before), well within the site limit of 120,000. Updated monitor/cancel hints.
- 🛠️**Updated** `gwaslab.process.array_for_submit.sh` — appends `--chrom ${SLURM_ARRAY_TASK_ID}` to the Python command when running as an array task; logs the chromosome in the job header.

## 2026-03-16 🆕 Reference file download utility
- 🆕**Added**: `gwaslab.download_refs.py` — utility script and complete inventory of all gwaslab reference files, using gwaslab's built-in `gl.download_ref()` function. Active entries (AFR, EAS, AMR, SAS for both hg19 and hg38) are downloaded; all other files already present at the reference directory are listed as comments and can be uncommented to (re-)download. Covered categories: 1KG population VCFs (all six populations, hg19 + hg38), HapMap3 EAF tables, 1KG SNPID→rsID conversion tables, dbSNP v151/v157 VCFs (very large, NCBI FTP), UCSC reference FASTA, recombination maps, and Ensembl/RefSeq GTF files. The `.tbi` index is fetched automatically alongside each VCF. Target directory defaults to `/path/to/references/gwaslab/`; override with `--ref-dir`. Prints a summary via `gl.check_downloaded_ref()` on completion.
- 🛠️**Updated**: `README.md` — added a `📥 Reference file management` section with a full reference-file inventory table (keyword, filename, default status), usage instructions, and a note about Dropbox/NCBI accessibility on HPC. Added `gwaslab.download_refs.py` to the HPC helper-scripts file table.

## 2026-03-16 🛠️ Updated environment.yml and installation instructions
- 🛠️**Updated**: `environment.yml` — overhauled to reflect the full dependency set required by `gwaslab.process.py`. Upgraded Python from `3.11` to `3.12`. Moved all Python packages to the `pip:` block with pinned or bounded versions: `numpy>=1.21.2,<2`, `adjusttext==0.8`, `matplotlib>=3.8,<3.9`, `pandas>=1.3,!=1.5`, `pysam==0.22.1`, `scikit-allel>=1.3.5`, `scipy>=1.12`, `seaborn>=0.12`, `h5py>=3.10.0`, `pyarrow`, `polars>=1.27.0`, `sumstats-liftover==1.1.0`, `jupyter==1.0.0`, `gwaslab`, `pyliftover`, `tqdm`. `bcftools` retained as a conda dependency (bioconda channel) rather than a pip package. Replaced `defaults` channel with `nodefaults` to avoid the Anaconda commercial repository, which is not permitted at many academic institutions; all packages are sourced exclusively from `conda-forge` and `bioconda`.
- 🛠️**Updated**: `README.md` — replaced the requirements and installation sections. Now documents Python 3.12 and `bcftools` as requirements; provides two installation paths (Option A: `mamba env create -f environment.yml`; Option B: manual `mamba create` + `pip install`); updated verification command to import `gwaslab` and `polars`; updated troubleshooting guidance for dependency conflicts and bioconda `bcftools`.

## 2026-03-16 🛠️ Two-tier resource model for gwaslab.process.submit_staged.sh
- 🛠️**Updated**: `gwas_list.txt` — added two new columns: `MEM_LIGHT` (COL10) and `TIME_LIGHT` (COL11) for the moderate pipeline stages (`process-normalize`, `process-check-ref`, `qc`). The existing `MEM` (COL8) and `TIME` (COL9) columns are unchanged and continue to control the heavy VCF-sweep stages (`process-infer-strand`, `process-assign-rsid`, `process-check-af`). Note added to header: `MEM_LIGHT` should be set higher for studies with many columns or complex allele structure (e.g. the AF multi-ancestry meta-analysis required 128G at `process-check-ref` despite having fewer variants than EUR studies that passed at 64G).
- 🛠️**Updated**: `gwaslab.process.submit_staged.sh` — replaced per-stage fixed defaults with two script-level fallback defaults (`MEM_LIGHT_DEFAULT=64G`, `MEM_HEAVY_DEFAULT=128G`). Per-study `MEM_LIGHT`/`TIME_LIGHT` are read from COL10/COL11 and applied to all light-tier stages (`process-normalize`, `process-check-ref`, `qc`); if absent the fallbacks are used. Report table now shows both tiers alongside the job chain.
- 🛠️**Updated**: Active entries in `gwas_list.txt` — `MEM_LIGHT`/`TIME_LIGHT` assigned per study: `32G/12h` for standard EUR studies; `64G/24h` for PAN and large EUR studies; `128G/24h` for AF (known to require higher memory at `process-check-ref`).

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

# CHANGES

This document tracks changes to the codebase. Each entry should include a brief description of the change, the files affected, and any relevant context or reasoning behind the change. This helps maintain a clear history of modifications and facilitates collaboration among developers.

## 2026-06-11 🐛 gwas_process.py — EAF lookup: normalise allele case before reference VCF matching (v1.4.32)
- **Root cause**: `check_and_fill_eaf()` built its `(pos, ref, alt)` lookup dict directly from the reference VCF (uppercase alleles: A/T/C/G) and then queried it using allele strings taken verbatim from the GWAS data. Older meta-analysis files (e.g. AholaOlli2017 cytokine GWAS) store alleles in lowercase (a/t/c/g). Python dict lookups are case-sensitive, so every query returned `None` — 0 out of ~9.9M EAF values were filled despite 9.8M rsIDs being present and the reference VCF chromosome names matching correctly (confirmed via `tabix -l`).
- **Diagnosis**: cross-referencing the `check_and_fill_eaf()` dict key construction (`fields[3]`, `fields[4]` from the VCF — always uppercase) against the per-row lookup (`str(row["ea"])`, `str(row["nea"])` — lowercase in affected files) confirmed a complete case mismatch. Chromosome-name format was ruled out as a secondary cause because the same reference VCF works correctly for other studies.
- **Fix**: two `.str.upper()` calls immediately after column renaming, before the per-chromosome loop — so all allele comparisons are uppercase-normalised regardless of how the source file encodes them.
- **Affected studies**: IL6_EUR_AholaOlli2017 (confirmed 0/9,901,590 EAF found); any study whose source file uses lowercase alleles will benefit from a rerun after this fix. Studies with uppercase alleles in source are unaffected.
- **Files**: `gwas_process.py` (v1.4.31 → v1.4.32).

## 2026-04-05 🐛 gwas_process.py — check_ref: prefer uncompressed FASTA to avoid pyfaidx OOM (v1.4.31)
- **Root cause**: `run_check_ref()` always used `hg{build}.fa.gz` (plain gzip) as the FASTA reference. pyfaidx cannot perform random-access on plain-gzip files — it decompresses and indexes the entire genome into RAM, which for hg38 (~3.2 billion bases as a Python in-memory structure) can exceed 100 GB regardless of variant count. Studies using hg38 as the target build (BUILD=38, or BUILD=19+liftover) hit this limit during `process-check-ref` even with only hundreds of thousands of per-chromosome variants.
- **Fix**: both the per-chromosome path (`run_check_ref()`) and the legacy whole-genome path now prefer `hg{build}.fa` (uncompressed) over `hg{build}.fa.gz`. With the uncompressed FASTA and its `.fai` index, pyfaidx uses true O(1) random access with memory proportional only to the variants being processed. pyfaidx creates the `.fai` automatically on first use if absent (one-time cost). `.fa.gz` is retained as a fallback with a warning.
- **Action required**: ensure `samtools faidx hg38.fa` has been run in REF_DIR so the `.fai` index exists (or let gwaslab/pyfaidx build it on first run). The uncompressed `hg38.fa` is already present in the reference directory alongside `hg38.fa.gz`.
- **Affected studies**: any study using the hg38 FASTA for check_ref, including all BUILD=19+liftover studies (IL6, CAD, HF, NICM, etc.) and BUILD=38 studies (CRP, Migraine, T1D).
- **Files**: `gwas_process.py` (v1.4.30 → v1.4.31).

## 2026-04-04 🐛 gwas_process.py — LDSC skips output when EAF is absent or all-NaN (v1.4.30)
- **Bug**: `write_ldsc()` applied `EAF > 0.01 & EAF < 0.99` even when EAF was all-NaN, silently removing every variant and writing a useless 0-variant LDSC file. Same root cause as the `apply_qc()` EAF/DAF bug (v1.4.29) — NaN comparisons always return False in pandas.
- **Different fix from apply_qc()**: for LDSC, EAF is genuinely required (it becomes the `Frq` column). Skipping the filter would produce an LDSC file with all-NaN frequencies, which is equally useless. Instead, `write_ldsc()` now detects all-NaN or absent EAF up-front, logs an ERROR explaining the cause and remediation (`--fill-eaf`), and returns without writing a file.
- **Studies affected**: `HF_EUR_Aragam2018` and `NICM_EUR_Aragam2018` used `--no-fill-eaf` and have no EAF in the source → will now clearly log the reason instead of writing an empty LDSC file. `T1D_EUR_Chiou2021` had wrong BUILD (19 instead of 38 for a GWAS Catalog harmonised file) → double liftover corrupted coordinates → EAF fill failed → near-0 LDSC variants; fix is BUILD=38 + full rerun.
- **Files**: `gwas_process.py` (v1.4.29 → v1.4.30).

## 2026-04-04 🐛 gwas_process.py — QC wipes all variants when EAF/DAF is all-NaN (v1.4.29)
- **Bug 1 (EAF)**: `build_qc_filter_expr()` always included the EAF filter regardless of whether EAF had any non-NaN values. In pandas, `NaN >= 0.005` evaluates to `False`, so an all-NaN EAF column caused every variant to fail the filter. All other filter terms (BETA, SE, INFO, DAF) were already guarded with `if col in cols else None`; EAF was not.
- **Bug 2 (DAF)**: `apply_qc()` passed `else 0` (not `else None`) when DAF column was absent, which would generate the impossible expression `DAF < 0 & DAF > 0` had the column not existed. Now passes `None`.
- **Root cause for CRP_EUR_Said2022**: EAF is mostly/entirely NaN in the harmonised input file and EAF fill could not recover it. The all-NaN EAF caused both the EAF filter and the DAF filter (DAF is derived from EAF) to remove all 10.6M variants.
- **Fix**: `apply_qc()` uses `_col_usable(col)` — column present AND has at least one non-NaN value — as the guard for every numeric filter (EAF, DAF, BETA, SE, INFO). `build_qc_filter_expr()` now accepts `eaf: float | None` and returns `None` when no usable filter criterion exists (all columns absent or all-NaN).
- **New**: when `expr` is `None` (no usable filters), all variants are retained with a WARNING. When QC removes 100% of variants, an ERROR is logged explicitly pointing to EAF/DAF as the likely culprit.
- **Files**: `gwas_process.py` (v1.4.28 → v1.4.29).

## 2026-04-04 ✨ gwas_process.py — OR + 95% CI columns in TSV output for case/control studies (v1.4.28)
- **New**: `reformat_output()` now detects case/control studies (N_cases present and non-zero) and appends three derived columns immediately after SE and P in the output TSV:
  - `OR` = exp(Beta)
  - `OR_lower_95CI` = exp(Beta − 1.96 × SE)
  - `OR_upper_95CI` = exp(Beta + 1.96 × SE)
- Because OR→BETA conversion happens at preprocess time (v1.4.27), `Beta` in all output files is already ln(OR), so these back-transformations are exact.
- Applies to both raw pre-QC (`.tsv.gz`) and QC-filtered (`.qc.tsv.gz`) outputs. COJO and LDSC outputs are unaffected (they require Beta in linear scale by design). Parquet/pickle internal formats are also unaffected.
- For quantitative traits (no N_cases) the output is unchanged.
- **Files**: `gwas_process.py` (v1.4.27 → v1.4.28).

## 2026-04-04 🐛 gwas_process.py — OR→BETA conversion at preprocess time (v1.4.27 supersedes v1.4.26)
- **Root cause identified**: `check_or_vs_beta()` (called during preprocess) only handled the *mislabelled* OR case (negative values → rename OR to BETA). When OR values were genuinely positive it logged "OR column looks valid" and left the column as `OR`. gwaslab preserves `OR` as-is through normalize → checkref → inferstrand → checkaf → merge, so `BETA` was never populated, and `write_cojo()` rightly skipped output.
- **Fix** (`check_or_vs_beta()`): genuine positive OR is now converted to `BETA = ln(OR)` and the `OR` column is dropped immediately at preprocess time. This means BETA flows correctly through the entire pipeline (gwaslab's `flip_allele_stats()` negates BETA, which is mathematically identical to taking ln(1/OR) = −ln(OR)), and COJO/LDSC/raw outputs all work without special casing.
- **Safety net retained** (`write_cojo()`): the OR→BETA fallback added in v1.4.26 remains in place as a guard for edge cases where OR survives to merge (e.g. old checkpoints written before this fix).
- **Affected studies**: any case/control study where the source file stores odds ratios (PGC iPSYCH ASD, PGC BIP OConnell2025, and likely ICH, AAA, TAA, IA, Migraine). Studies with `.withBETA_N` files (BIP_EUR_Stahl2019, CD/IBD/UC Liu2015) and UKBB Neale files (which publish log-OR as BETA) are unaffected.
- **Files**: `gwas_process.py` (v1.4.26 → v1.4.27).

## 2026-04-04 🐛 gwas_process.py — COJO OR→BETA conversion for case/control studies (v1.4.26)
- **Bug**: `write_cojo()` required a `BETA` column and silently skipped COJO output when the study used odds ratios (OR). Affected OR-based case/control studies (e.g. ASD_EUR_Grove2019, BIP_EUR) which had no `BETA` column after gwaslab processing.
- **Fix**: if `BETA` is absent but `OR` is present, `write_cojo()` now derives `BETA = ln(OR)` on a working copy of the DataFrame before building the COJO table. Variants with `OR ≤ 0` or non-finite OR are set to `NaN` and filtered out with a warning.
- **Safety**: a non-finite filter (`np.isfinite(BETA) & np.isfinite(SE)`) is applied to the COJO output regardless of whether OR→BETA conversion was performed, guarding against any upstream data issues.
- **Files**: `gwas_process.py` (v1.4.25 → v1.4.26).

## 2026-04-04 🐛 gwas_process.py — fix infer_ancestry hg38 reference not found (v1.4.25)
- **Bug**: `run_infer_ancestry()` used gwaslab's key-based file lookup (`1kg_hm3_hg38_eaf`) without first telling gwaslab where to look. gwaslab searches its own default cache; if the file was placed in `REF_DIR` manually (or via our download helper rather than gwaslab's internal helper), the lookup fails with `Reference file '1kg_hm3_hg38_eaf' not found` even when the file is physically present.
- **Fix**: `run_infer_ancestry()` now accepts an optional `ref_dir` keyword argument. When provided and the directory exists, `gwaslab.bd.bd_download.set_default_directory(ref_dir)` is called before `infer_ancestry()` so gwaslab resolves the HapMap3 EAF reference (`PAN.hapmap3.hg38.EAF.tsv.gz` / `PAN.hapmap3.hg19.EAF.tsv.gz`) from the correct location. Both call sites pass `ref_dir=args.ref`.
- **Files**: `gwas_process.py` (v1.4.24 → v1.4.25).

## 2026-04-04 🐛 make_chrpos_hdf5.sh — fix gwas2cojo.conf not found in SLURM spool
- **Bug**: SLURM copies the script to its spool directory before execution, making `BASH_SOURCE[0]` resolve to `/var/spool/slurmd/jobN/slurm_script` rather than the original script location. The conf lookup `${SCRIPT_DIR}/../gwas2cojo.conf` therefore failed with `ERROR: gwas2cojo.conf not found`.
- **Fix**: conf resolution now follows a three-step fallback:
  1. `GWAS2COJO_CONF` env var — explicit override
  2. `${SLURM_SUBMIT_DIR}/gwas2cojo.conf` — SLURM always exports `SLURM_SUBMIT_DIR` as the directory where `sbatch` was called; submitting from the gwas2cojo root just works
  3. `BASH_SOURCE` relative path — fallback for direct local invocation
- **New**: after sourcing the conf, `ROOTDIR=$(dirname "${PYTHON_SCRIPT}")` derives the installation root from the already-known `PYTHON_SCRIPT` path, so the Python helper is always called as `${ROOTDIR}/utility_scripts/make_chrpos_hdf5.py` regardless of spool location.
- **File**: `utility_scripts/make_chrpos_hdf5.sh`.

## 2026-04-04 🐛 check.py — array stage chromosome count uses split ground truth (v1.2.7)
- **Fix**: array stages (`checkref`, `inferstrand`, `assignrsid`, `checkaf`) now use the chromosome count reported by the split stage (`n_split_chr`) as the expected total (`eff_total`) rather than counting SLURM log files. SLURM always arrays over all 26 tasks regardless of how many chromosomes are in the data; tasks for absent chromosomes exit 0 without a "done" marker, previously causing false `⚠ 23/26` warnings for datasets with only autosomes + chrX.
- **Behaviour**: if split reports N chromosomes and all N array tasks complete, status shows `✓ done` with metric `N chromosomes complete`. A warning is only raised when `n_done < n_split_chr` (genuine missing or failed tasks).
- **Files**: `gwas_process.check.py` (v1.2.6 → v1.2.7).

## 2026-04-03 ✨ gwas_process.py — --add-chrpos flag + make_chrpos_hdf5.py utility (v1.4.24)
- **New** (`gwas_process.py`): `--add-chrpos` flag assigns CHR and POS from rsID at the preprocess stage for datasets that contain only rsIDs (e.g. MVP CAD files). Must be added per-study via `EXTRA_FLAGS` (COL12) in `gwas_list.txt`; not set globally.
- **New** (`gwas_process.py`): `assign_chrpos_from_hdf5()` implements the lookup directly on the pandas DataFrame using pre-built per-chromosome HDF5 files in `--ref`. Uses the same modulo-10 group structure as gwaslab's `rsid_to_chrpos2()`. Parallel lookup via `ThreadPoolExecutor`. When CHR is absent, searches all chromosome files; when CHR is present, restricts to matching files.
- **New** (`utility_scripts/make_chrpos_hdf5.py`): one-time setup script wrapping `gl.process_vcf_to_hfd5()`. Reads `REF_DIR` from `gwas2cojo.conf`, discovers the best available dbSNP VCF (`v157` preferred, `v151` fallback) via gwaslab's `get_path()`, and writes HDF5 files into `REF_DIR`. Options: `--ref-dir`, `--build hg19|hg38|all`, `--threads`, `--complevel`, `--overwrite`.
- **Usage**: run `make_chrpos_hdf5.py --build hg19` once, then add `--add-chrpos` to `EXTRA_FLAGS` for affected studies. Multiple flags in `EXTRA_FLAGS` are space-separated (e.g. `--add-chrpos --keep-multiallelic`).
- **Files**: `gwas_process.py` (v1.4.23 → v1.4.24), `utility_scripts/make_chrpos_hdf5.py` (new, v1.0.0).

## 2026-04-03 ⚡ gwas_process.py — vectorise EAF fill + check.py EAF reporting (v1.4.23 / check v1.2.6)
- **Fix** (`gwas_process.py`): `check_and_fill_eaf()` rewritten to fetch per chromosome via tabix rather than issuing one tabix query per variant. For a 3.5 M-variant file across 22 chromosomes this reduces I/O from ~3.5 M individual tabix calls to ~22, cutting runtime from hours to seconds and eliminating preprocess TIMELIMIT kills.
- **How**: for each chromosome, one tabix range-fetch covers all positions in that contig; results are loaded into a `(pos, ref, alt) → AF` dict; per-variant AF is resolved by dict lookup (O(1)) with automatic allele-flip when effect/other alleles are swapped.
- **New** (`gwas_process.check.py`): `_metrics_preprocess()` now parses EAF fill log lines from the preprocess `.out` file and surfaces them in the check summary row:
  - `EAF: complete` — EAF column present and fully populated
  - `EAF: not filled` — `--fill-eaf` not set or suppressed
  - `EAF filled N from ref (M still missing)` — partial or full fill from reference VCF
- **Files**: `gwas_process.py` (v1.4.22 → v1.4.23), `gwas_process.check.py` (v1.2.5 → v1.2.6).

## 2026-03-27 🔁 rename: gwaslab.process.* → gwas_process.*
- **Rename**: all pipeline scripts renamed from `gwaslab.process.<name>` to `gwas_process.<name>` for consistency and brevity:
  - `gwas_process.py` → `gwas_process.py`
  - `gwas_process.check.py` → `gwas_process.check.py`
  - `gwas_process.cleanup.sh` → `gwas_process.cleanup.sh`
  - `gwas_process.submit.sh` → `gwas_process.submit.sh`
  - `gwas_process.submit_staged.sh` → `gwas_process.submit_staged.sh`
  - `gwas_process.array_for_submit.sh` → `gwas_process.array_for_submit.sh`
- **Updated**: all internal cross-references, `VERSION_NAME`, `prog=`, log file suffix (`.gwaslab_process.log` → `.gwas_process.log`), and submit-log filename prefix updated accordingly.
- **Files**: all six scripts above + `CHANGES.md`.

## 2026-03-26 🔧 check.py — outputs row + ancestry display fix (check v1.2.4)
- **Fix**: `_parse_ancestry_check()` now correctly handles `Match: unknown  ⚠ SKIPPED` log lines (emitted when EAF is absent/all-NaN and `infer_ancestry` is skipped). Previously the `UNKNOWN` match value fell through to `status: "unknown"`, showing `ancestry: unknown (status unknown)`. Now mapped to `status: "skipped"` and displayed as `ancestry: not inferred — EAF unavailable (provided=POP)`.
- **New**: `_metrics_outputs(text)` extracts COJO and LDSC output variant counts from the merge stage log (`[SAVE] COJO → ...` and `[SAVE] LDSC → ...` patterns).
- **New**: A `└ outputs` sub-row is printed directly after the merge row when either COJO or LDSC (or both) outputs were written, showing variant counts (e.g. `COJO 6,912,451  |  LDSC 1,103,847`). COJO count removed from the merge row itself.
- **Files**: `gwas_process.check.py` (v1.2.3 → v1.2.4).

## 2026-03-25 ✨ LDSC-ready output via --ldsc flag (v1.4.21)
- **New**: `write_ldsc()` function produces an LDSC-ready munged summary statistics file from the QC-filtered data. Applies the standard LDSC pre-filtering pipeline on an isolated deep copy (original QC object unchanged):
  1. `filter_hapmap3()` — HapMap3 variants only
  2. `filter_palindromic(mode="out")` — all A/T and C/G SNPs removed (LDSC cannot handle strand ambiguity)
  3. `exclude_hla()` — HLA region excluded (chr6:25–34 Mb)
  4. `filter_region_out(high_ld=True, build=output_build)` — other high-LD regions excluded
  5. `filter_value('INFO > 0.9 & EAF > 0.01 & EAF < 0.99')` — quality thresholds (INFO filter skipped if column absent)
  6. `to_format(fmt="ldsc")` — gwaslab ldsc format: SNP (rsID), A1, A2, Beta/OR, Frq, INFO, N, P, Z, CHR, POS
- **New**: `--ldsc` argparse flag (analogous to `--cojo`); enabled by default in both submission scripts.
- **Output**: `{stem}.qc.ldsc.tsv.gz` alongside the existing `.qc.tsv.gz` and `.cojo.gz`.
- **Robustness**: each filter step wrapped in try/except — a missing reference file or failed step skips that step and logs a warning without aborting the pipeline.
- **Files**: `gwas_process.py` (v1.4.20 → v1.4.21), `gwas_process.array_for_submit.sh`, `gwas_process.submit_staged.sh`.

## 2026-03-25 ✨ ancestry inference check at QC stage (v1.4.20 / check v1.2.3)
- **New** (`gwas_process.py`): `run_infer_ancestry()` calls `gwas_obj.infer_ancestry()` on the QC-filtered data, comparing the declared `--population` against the Fst-inferred super-population from the HapMap3 pan-ancestry EAF reference (`1kg_hm3_hg19/hg38_eaf`). Run at the end of the QC block in both the `merge` stage and the `--stage all` path.
- **New** (`gwas_process.py`): `--no-infer-ancestry` flag skips the ancestry check (enabled by default). Logged under Toggles as `infer_ancestry=True/False`.
- **Logging**: emits a canonical `[ANCESTRY CHECK] Provided: X | Inferred: Y | Match: True/FALSE` line (WARNING level on mismatch) parseable by the check script.
- **Output**: result saved to `{stem}.ancestry_check.json` in the output directory for archival.
- **New** (`gwas_process.check.py` v1.2.3): `_parse_ancestry_check()` parses the `[ANCESTRY CHECK]` log line from the merge or qc stage output. Result displayed in the study header line. Mismatches shown as `⚠ MISMATCH` in the header and `⚠ ANCESTRY MISMATCH — re-check population label` in the overall summary. `--errors-only` also surfaces ancestry mismatches.
- **Files**: `gwas_process.py` (v1.4.19 → v1.4.20), `gwas_process.check.py` (v1.2.2 → v1.2.3).

## 2026-03-25 ✨ gwaslab.process.py — comprehensive STATUS filter + --filter-palindromic (v1.4.19)
- **New**: `_apply_status_filter()` helper replaces the previous single digit_7 check with a comprehensive STATUS-based filter covering all problematic flag classes:
  - Build prefix 97/98: `UnknownGenome` / `UnmappedVariant` (e.g. liftover failures)
  - Digit 4 in [5,6,7,8]: CHR or POS invalid/unknown (safety net; most handled by `basic_check(remove=True)`)
  - Digit 5 in [5,6,7]: allele indistinguishable, invalid notation, or unknown
  - Digit 6 = 8: not on reference genome (safety net; most removed by `check_ref` internally)
  - Digit 7 in [7,8]: `infer_strand2` indistinguishable (7) or no match/no info (8) — previously only 8 was filtered
- **New**: `--filter-palindromic` flag calls `filter_palindromic(mode="out")` to remove ALL A/T and C/G SNPs at QC. Disabled by default — the STATUS filter is more precise (resolved palindromics at asymmetric MAF are retained; only unresolvable ones removed via digit_7 [7,8]). Use for strict meta-analysis strand-safety.
- **Logging**: STATUS filter reports per-class counts so the breakdown is visible in the log.
- **Note**: Digit 3 (SNPID/rsID format) issues are intentionally not filtered — they represent ID format problems only; CHR:POS and alleles are still valid and usable.
- **Files**: `gwas_process.py` (v1.4.18 → v1.4.19).

## 2026-03-25 ✨ gwaslab.process.py — STATUS-based filter at QC stage (v1.4.18)
- **New**: `apply_qc()` now runs a STATUS digit-7 filter after the numeric threshold pass. Variants where `infer_strand2` could not resolve the strand (STATUS digit_7 == 8 — palindromic SNPs at MAF~0.5, or indel allele mismatches) are removed before saving QC output.
- **Background**: gwaslab's `check_ref` already internally removes variants with digit_6 == 8 (allele absent from FASTA reference). `check_af2` does not use STATUS — it populates the DAF column, which is covered by `--daf-max`. The only STATUS flag that survives to output without removal is digit_7 == 8 from `infer_strand2`. No `filter_status()` method exists in this gwaslab version; the filter is implemented directly via integer arithmetic (`STATUS % 10 == 8`).
- **Logging**: separate counts for numeric filter and STATUS filter; total after all QC filters logged at end.
- **Files**: `gwas_process.py` (v1.4.17 → v1.4.18).

## 2026-03-25 ✨ gwaslab.process.py — add normalize_allele + basic_check(remove=True) (v1.4.17)
- **New**: `normalize_allele(threads=n_cores)` inserted between `basic_check()` and `remove_dup()` in both `run_normalize()` and the `--stage all` path. This standardises indel notation (trim shared prefix/suffix, uppercase, left-align) before deduplication so that variants expressed differently across studies but representing the same position are correctly identified as duplicates.
- **Change**: `basic_check()` now called with `remove=True` (previously no arguments), so variants with invalid chromosome codes, positions, or allele strings are dropped at source rather than propagating through the pipeline.
- **Files**: `gwas_process.py` (v1.4.16 → v1.4.17).

## 2026-03-25 ✨ gwaslab.download_refs.py — add chromosome X reference support (v1.2.0)
- **New**: Added `_ANCESTRY_X_VCFS` dict covering all 6 ancestries × 2 builds for the 1KG chrX VCFs (`1kg_{eur,pan,afr,eas,amr,sas}_x_hg19/hg38`).
- **New**: Added `1kg_dbsnp151_hg19_x` and `1kg_dbsnp151_hg38_x` chrX SNPID→rsID conversion tables to `_BUILD_FILES` (alongside existing autosomal `_auto` tables).
- **New**: `--no-x` flag to skip all chrX downloads (population VCFs + rsID tables); default behaviour is to include chrX.
- **New**: `build_download_list()` gains `include_x` parameter; `main()` prints `Include chrX` status line.
- **Files**: `gwaslab.download_refs.py` (v1.1.0 → v1.2.0).

## 2026-03-24 🐛 cleanup.sh — fix --config rejecting process substitution
- **Bug**: `[[ ! -f "${CONFIG_FILE}" ]]` uses `-f` which only matches regular files; process substitution (`<(...)`) passes a named pipe (`/dev/fd/N`) which fails the test, producing `ERROR: config file not found: /dev/fd/63`.
- **Fix**: changed to `[[ ! -e "${CONFIG_FILE}" ]]` (`-e` matches any file type including pipes), so `--config <(grep ...)` now works as expected.
- **File**: `gwas_process.cleanup.sh`

## 2026-03-23 🐛 gwaslab.process.py — alias table corrections and gwas_list.txt fixes (v1.4.16)
- **Fix**: added `"freq_a"` to EAF aliases in `SUMSTATS_ALIASES` — covers TAG consortium files (`tag.*.tbl.withN.txt.gz`) which use `FRQ_A` for effect-allele frequency. Previously EAF was always NaN for these studies and had to be filled from the reference VCF.
- **Fix**: added `"imp_qual"` to INFO aliases — covers MVP PAD file (`CLEANED.MVP.EUR.PAD.results.anno.nodup.txt.gz`) which stores imputation quality as `IMP_QUAL`.
- **Fix**: added `"rs_id"`, `"dbsnp_rs_id"`, `"dbsnp_id"` to rsID aliases — broadens coverage for MVP and dbSNP-derived header variants.
- **Fix** (`gwas_list.txt`): `MI_PAN_withmultiallelic` path was missing the leading `/` — would have failed at file open.
- **Fix** (`gwas_list.txt`): `AF_PAN` build was `19` but the TOPMed Freeze 5 file uses `position_b38` (hg38 coordinates) — corrected to `build=38` to prevent double-liftover.
- **Fix** (`gwas_list.txt`): `MEM_LIGHT` was `64GB` for `PAD_EUR_MVP`, `PAD_EUR_FINNGEN`, `PAD_EUR_UKB` — SLURM requires `64G`; corrected.
- **Change** (`gwas_list.txt`): TAG study names standardised to `TAG_EUR_*` naming convention for consistency with other phenotype–ancestry naming in the list.
- **Files**: `gwas_process.py` (v1.4.15 → v1.4.16), `gwas_list.txt`.

## 2026-03-20 ✨ Add --no-fill-eaf per-study override flag (v1.4.15)
- **Problem**: The submit script passes `--fill-eaf` globally for all studies. Studies with no EAF column trigger a per-variant tabix lookup across the full VCF for every variant (O(n)), which is prohibitively slow for large files (e.g. 7.7M variants × tabix = many hours).
- **Fix**: Added `--no-fill-eaf` flag that suppresses the EAF lookup even when `--fill-eaf` is present. Intended for use as a per-study `EXTRA_FLAGS` override in `gwas_list.txt`.
- **Usage**: Add `;--no-fill-eaf` as COL12 in `gwas_list.txt` for the affected study. EAF will be filled properly at the `process-check-af` stage from the 1KG VCF anyway.
- **Logging**: `--no-fill-eaf` is logged in the Toggles line. A separate info message confirms suppression when both flags are present.
- **Files**: `gwas_process.py` v1.4.15.

## 2026-03-20 🐛 Header typo: missing CHR/NEA aliases and KeyError in remove_dup (v1.4.14)
- **Root cause**: The Suzuki2024 T2DGGI file uses `Chromsome` (typo, missing 'o') for chromosome and `NonEffectAllele` (no underscore) for the non-effect allele. Neither matched existing aliases, so gwaslab never received a CHR or NEA column.
- **Symptom 1**: `KeyError: Index(['CHR'])` in `run_normalize()` at the `duplicated(subset=["CHR","POS"])` multi-allelic count — crashed before `remove_dup` was even called.
- **Fix 1**: Added `"chromsome"` to CHR aliases in `SUMSTATS_ALIASES` (typo-tolerant match).
- **Fix 2**: Added `"noneffectallele"` (no underscore) to NEA aliases in `SUMSTATS_ALIASES`.
- **Fix 3**: Guarded the `duplicated(subset=["CHR","POS"])` call in both `run_normalize()` and `run_processing()` with a column-existence check (`_has_chr_pos`) so a missing CHR column logs 0 multi-allelics rather than raising `KeyError`.
- **Note**: No SNPID column in this file is not a blocker — gwaslab derives CHR:POS:NEA:EA IDs via `fix_id` once CHR and NEA are correctly mapped.
- **Files**: `gwas_process.py` v1.4.14.

## 2026-03-19 ✨ gwaslab.download_refs.py — conf-file integration, --build and --ancestry arguments
- **Feature**: reference directory now defaults to `REF_DIR` from `gwas2cojo.conf` (parsed next to the script) instead of a hardcoded placeholder. A warning is printed if the conf is absent or `REF_DIR` is unset.
- **Feature**: new `--build` argument (`all` / `hg19` / `hg38`, default: `all`) filters downloads to the requested genome build(s).
- **Feature**: new `--ancestry` argument (`EUR` / `PAN` / `AFR` / `EAS` / `AMR` / `SAS` / `all`, default: `EUR`) selects which 1KG population VCF(s) to download. Default is `EUR` to avoid accidentally triggering all 12 large VCF downloads.
- **Change**: the flat `TO_DOWNLOAD` list replaced by structured `_ANCESTRY_VCFS` and `_BUILD_FILES` dictionaries; `build_download_list()` assembles the final keyword list at runtime based on the selected builds and ancestries. Non-ancestry-specific files (dbSNP, FASTA, recombination maps, GTFs, HapMap3 EAF, SNPID→rsID tables) are always included for the requested build(s).
- **Files**: `gwaslab.download_refs.py`, `README.md`.

## 2026-03-19 ✨ gwaslab.process.py — detect OR column mislabelled as BETA and auto-rename (v1.4.13)
- **Feature**: new `check_or_vs_beta()` function called during preprocess (after `standardise_columns`). If the standardised `OR` column contains any negative values it cannot be a true odds ratio — the source file has mislabelled a BETA/log-odds column as `OR`. The function renames `OR` → `BETA` with a warning log line showing the count and percentage of negative values, and processing continues normally. If both `OR` and `BETA` are already present the check is skipped.
- **Example**: `tag.logonset.tbl.withN.txt.gz` has a column named `OR` containing effect sizes like `-0.0054`, `-0.0049`, which are clearly log-odds / BETA values. Without this fix the mislabelled column passed through to `run_check_ref` where `flip_allele_stats` tried to compute `1 / OR` and either hit `FloatingPointError` (OR = 0) or silently produced nonsensical results.
- **Files**: `gwas_process.py` (v1.4.12 → v1.4.13).

## 2026-03-19 🐛 gwaslab.process.py — FloatingPointError in flip_allele_stats when OR = 0 (v1.4.12)
- **Bug**: `run_check_ref` crashed with `FloatingPointError: divide by zero encountered in divide` inside gwaslab's `flip_by_inverse` when flipping OR-based studies (e.g. TAG_LogOnset). gwaslab computes `OR = 1 / OR` for flipped variants; if any OR value is 0 (missing data stored as zero rather than NaN), this raises a `FloatingPointError`. The error affected all 22 chromosomes (44 total errors = 2 per chromosome).
- 🐛**Fixed** `gwas_process.py` — `run_check_ref` now checks for an `OR` column before calling `flip_allele_stats`. Any rows where `OR = 0` are dropped with a warning log line before the flip, preventing the divide-by-zero. OR = 0 is not a biologically valid value; these are treated as missing data.
- **Files**: `gwas_process.py` (v1.4.11 → v1.4.12).

## 2026-03-19 🐛 gwaslab.process.py — plot_mqq crashes with TypeError when EAF is entirely missing (v1.4.11)
- **Bug**: `run_merge` → `plot_full_dataset` crashed with `TypeError: cannot unpack non-iterable NoneType object` when the dataset had no valid EAF values (e.g. DIAMANTE-TA PAN file has no EAF column). gwaslab's `_mqqplot` returns `None` instead of `(plot, log)` when it finds no plottable data, and the call site did not guard against this. The `plot_daf` call already had a `try/except`, but `plot_mqq` (both pre-QC and QC) did not.
- 🐛**Fixed** `gwas_process.py` — both `plot_mqq` loops (pre-QC in `plot_full_dataset` and QC in `plot_qc_dataset`) are now wrapped in `try/except TypeError` that logs a warning and skips the plot rather than crashing the pipeline. All other outputs (parquet, TSV.GZ, COJO, leads) are unaffected.
- **Files**: `gwas_process.py` (v1.4.10 → v1.4.11).

## 2026-03-19 🐛 gwaslab.process.check.py — normalize "after dedup" count not shown for default mode (v1.2.2)
- **Bug**: the normalize metric showed `liftover → 38` only (no variant count) for studies run with the default `mode="md"`. Root cause: `gwas_process.py` v1.4.9 changed the log message from `"After duplicate removal"` to `"After multi-allelic and duplicate variant removal"`, but the regex in `_metrics_normalize` still matched only the old wording. Studies run with `--keep-multiallelic` (`mode="d"`) still wrote the old message, so they showed the count while default-mode studies did not.
- 🐛**Fixed** `gwas_process.check.py` — `_metrics_normalize()` regex broadened to `After (?:multi-allelic and )?duplicate(?:\s+variant)? removal:` to match both message variants.
- **Files**: `gwas_process.check.py` (v1.2.1 → v1.2.2).

## 2026-03-19 🐛 gwaslab.process.check.py — "after dedup" variant count missing thousands separator (v1.2.1)
- **Bug**: the normalize metric displayed the post-dedup variant count without comma separators (e.g. `20073068 after dedup`) because `_first()` returns the raw matched string and the log line itself omits commas.
- 🐛**Fixed** `gwas_process.check.py` — `_metrics_normalize()` now converts the matched string to `int` and reformats it with `{:,}` before appending to the metric string, yielding `20,073,068 after dedup`.
- **Files**: `gwas_process.check.py` (v1.2.0 → v1.2.1).

## 2026-03-19 🐛 gwaslab.process.check.py — false ⚠ 22/26 chr when submit script arrays over 26 but split only produced 22 (v1.2.0)
- **Bug**: array stages (checkref, inferstrand, assignrsid, checkaf) reported `⚠ 22/26 chr` and set `any_error = True` for autosome-only datasets. The submit script always arrays over all 26 chromosomes; for non-autosomal chromosomes the job finds no input data and finishes without writing a `[SAVE]` marker, so `_is_done()` returned `False` for those 4 jobs. The code then saw `n_done=22, n_total=26` and flagged a warning even though every autosomal chromosome completed successfully.
- 🐛**Fixed** `gwas_process.check.py`:
  - After processing the `split` stage, the chromosome count is now stored in `n_split_chr`.
  - A new `split_autosome_only` flag is set when `n_split_chr == 22`; it is OR-ed with the existing `n_non_auto == 0` check into a combined `autosome_only` flag.
  - Array stages now track `n_auto_done` (autosomal chromosomes that completed) separately from `n_done` (all chromosomes). When `autosome_only`, the effective counts `eff_done / eff_total` are `n_auto_done / 22`, so non-autosomal "not done" files are ignored.
  - `checkref` aggregation is restricted to autosomal chromosome texts when `split_autosome_only` and non-autosomal log files exist, preventing empty files from skewing match-rate stats.
- **Files**: `gwas_process.check.py` (v1.1.0 → v1.2.0).

## 2026-03-19 ✨ gwaslab.process.check.py — wider metric column, full "unmatched" display, and autosome-only detection (v1.1.0)
- **Fix**: metric column widened from 46 to 58 characters (table width 100 → 112) so the full `unmatched N,NNN,NNN` value is no longer truncated to `unmat…` in the checkref row.
- **Feature**: array stages (checkref, inferstrand, assignrsid, checkaf) now distinguish between truly incomplete runs and datasets that contain only the 22 autosomes. When all 22 autosomal chromosomes completed and no non-autosomal files exist, the status is `✓ done` instead of the misleading `⚠ 22/26 chr`, and `any_error` is no longer set — making real failures much easier to spot.
- **Feature**: inferstrand / assignrsid / checkaf metric in autosome-only mode shows `22 autosomes complete` instead of `22/22 complete`.
- **Feature**: split metric now appends `, no non-autosomal` when exactly 22 chromosomes are present.
- **Files**: `gwas_process.check.py` (v1.0.0 → v1.1.0).

## 2026-03-19 ✨ Per-study extra flags via COL12 in gwas_list.txt
- **Feature**: `gwas_list.txt` now supports an optional 12th semicolon-delimited field (`EXTRA_FLAGS`) for per-study flags passed verbatim to `gwas_process.py`. Use `.` as a no-op placeholder. Multiple flags are space-separated within the field.
- **Example**: append `;--keep-multiallelic` to a study line to retain multi-allelic variants for that study only, while all other studies use the default `mode="md"` removal.
- **Example**: `;--keep-multiallelic --no-figures` to combine multiple flags.
- **Files**: `gwas_process.array_for_submit.sh` (reads COL12, appends to CMD array); `gwas_process.submit_staged.sh` (parses COL12 and documents it; the field passes through to the worker via `LINE`).
- **Logging**: the worker script echoes `Extra flags : <value>` in the job header for traceability.

## 2026-03-19 ✨ Extended column alias coverage for three new GWAS header formats (gwaslab.process.py v1.4.10)
- **Feature**: Added aliases for three additional GWAS summary statistics header formats:
  - **Format 1** (`Tested_Allele` / `Freq_Tested_Allele_in_HRS`): `tested_allele` → EA; `freq_tested_allele_in_hrs` → EAF.
  - **Format 2** (meta-analysis fixed-effects): `chromosome(b37)` → CHR; `position(b37)` → POS; `chrposid` → SNPID; `fixed-effects_beta` → BETA; `fixed-effects_se` → SE; `fixed-effects_p-value` → P.
  - **Format 3** (GWAS Catalog harmonised `hm_*` columns): `hm_variant_id` / `variant_id` → SNPID; `hm_rsid` → rsID; `hm_chrom` → CHR; `hm_pos` → POS; `hm_effect_allele` → EA; `hm_other_allele` → NEA; `hm_beta` → BETA; `hm_effect_allele_frequency` → EAF.
- **Design**: `hm_*` aliases are placed **before** their bare equivalents in each list so that when a harmonised GWAS Catalog file contains both `hm_effect_allele` and `effect_allele`, the harmonised column is preferred by `resolve_column()`.

## 2026-03-19 ✨ --keep-multiallelic flag and multi-allelic count logging (gwaslab.process.py v1.4.9)
- **Feature**: new `--keep-multiallelic` flag. By default `remove_dup` runs with `mode="md"` (remove duplicates **and** multi-allelic variants). With `--keep-multiallelic` it runs with `mode="d"` (duplicates only), leaving multi-allelic sites in the dataset. Useful when the GWAS reports genuine multi-allelic signals or for exploratory analysis before committing to a COJO run.
- **Feature**: before calling `remove_dup`, the number of variants at multi-allelic positions (same CHR:POS, different alleles) is now counted via `duplicated(subset=["CHR","POS"], keep=False)` and included in the post-removal log line: `After multi-allelic and duplicate variant removal: N variants remain (X removed; Y variants were at multi-allelic positions).`
- **Change**: log message changed from `"After duplicate removal"` to `"After multi-allelic and duplicate variant removal"` to accurately describe what was removed.
- 🛠️**Updated**: `run_normalize()` and `run_processing()` both receive the new `keep_multiallelic` kwarg; both call sites in `main()` pass `keep_multiallelic=args.keep_multiallelic`.

## 2026-03-19 ✨ gwaslab.process.cleanup.sh — granular pickle retention flags; remove qc.pkl by default (cleanup.sh)
- **Change**: `KEEP_QC_PKL` default flipped from `1` → `0`. The `*.qc.pkl` contains only `self.data` (identical to `.qc.parquet`), `self.log` (redundant with archived SLURM logs), and gwaslab internal state flags — nothing needed for downstream analysis. Pass `--keep-qc-pkl` to retain it.
- **Feature**: `*.normalize.pkl` moved out of the unconditional removal list and into its own conditional block, controlled by `--keep-normalize-pkl` (default: remove). Useful if you want to reload the normalized Sumstats object without re-running preprocess + normalize.
- **Change**: `--remove-qc-pkl` flag removed (now the default); replaced by `--keep-qc-pkl` to opt in to retention.
- **Change**: raw pkl loop now also skips `*.normalize.pkl` (handled by its own block) in addition to `*.qc.pkl`.
- 🛠️**Updated**: header comment and usage examples updated accordingly.

## 2026-03-19 🐛 gwaslab.process.cleanup.sh — per-chromosome intermediate parquets not removed (cleanup.sh)
- **Bug**: per-chromosome intermediate parquets (`*.chr*.normalize.parquet`, `*.chr*.checkref.parquet`, `*.chr*.inferstrand.parquet`, `*.chr*.assignrsid.parquet`, `*.chr*.checkaf.parquet`) were not included in the cleanup patterns. These files are the stage-to-stage handoffs written by `save_chrom_parquet()` for each of the 26 chromosome array tasks, and collectively represent the largest share of intermediate disk usage (e.g. 26 × 5 stages × ~28 MB = ~3.6 GB per study for a large GWAS).
- 🐛**Fixed** `gwas_process.cleanup.sh` — added all five `*.chr*.{stage}.parquet` glob patterns to the `patterns` array in `cleanup_study()`. Updated the header comment to document them.

## 2026-03-19 ✨ gwaslab.process.cleanup.sh — archive SLURM log files into study/logs/ (cleanup.sh)
- **Feature**: after removing intermediate checkpoints, the cleanup script now moves all SLURM `*.out` / `*.err` files belonging to the study from `LOG_DIR` (default: `OUT_BASE`) into `${OUT_BASE}/<STUDY>/logs/`. This keeps the submit directory tidy and preserves logs in a study-specific location for future use with `gwas_process.check.py`.
- New flags: `--no-archive-logs` (skip archiving), `--log-dir PATH` (override source directory when SLURM logs land elsewhere).
- After archiving, the script prints the exact `gwas_process.check.py` command to inspect the archived logs.
- Log archiving is enabled by default; `--dry-run` mode previews what would be moved without touching files.

## 2026-03-19 🐛 gwaslab.process.cleanup.sh — wrong output directory path and fragile glob (cleanup.sh)
- **Bug**: `gwas_process.cleanup.sh` was silently doing nothing in all three modes (`--study`, `--all`, `--config`). Root cause: all three study-directory paths were constructed as `${OUT_BASE}/${STUDY_NAME}/GWASCatalog`, but `/GWASCatalog` is only appended by `gwas_process.py` when `--output` is *not* passed on the command line. The pipeline always passes `--output "${OUT_BASE}/${GWAS_NAME}"` explicitly (via `array_for_submit.sh`), so `output_loc = args.output` — no `/GWASCatalog` suffix. Every directory existence check therefore failed and cleanup was skipped without any error.
- 🐛**Fixed** `gwas_process.cleanup.sh`:
  - `--study` mode: `${OUT_BASE}/${STUDY_NAME}/GWASCatalog` → `${OUT_BASE}/${STUDY_NAME}`
  - `--all` mode: glob `"${OUT_BASE}"/*/GWASCatalog` → `"${OUT_BASE}"/*/`; `basename "$(dirname ...)"` → `basename "${study_dir}"`
  - `--config` mode: `${OUT_BASE}/${GWAS_NAME}/GWASCatalog` → `${OUT_BASE}/${GWAS_NAME}`
  - Fragile unquoted glob `local raw_pkl_pattern="${study_dir}/"*.pkl` (expands at assignment time) → quoted `"${study_dir}/*.pkl"` (expands at `for` loop time)

## 2026-03-19 🐛 gwaslab.process.check.py — AttributeError on optional regex group (check.py v1.0.1)
- **Bug**: `AttributeError: 'NoneType' object has no attribute 'strip'` in `_first()` when called with a regex containing an optional capturing group (`(...)?`). The outer `re.search()` matched (so `m` was not `None`), but `m.group(1)` was `None` because the optional group did not participate in the match. The `m.group(group).strip() if m else default` guard only checked for a missing match, not for a `None` group value.
- 🐛**Fixed** `gwas_process.check.py` — `_first()` now checks `val = m.group(group)` separately and returns `default` if `val is None`. The broken `\[SAVE\] QC Parquet` regex with an optional group was also removed (it was dead code — the result was never used; `qc_n` via the `After QC` pattern was the operative extraction). Bumped to v1.0.1.

## 2026-03-18 🧰 Added gwaslab.process.check.py — pipeline run-status checker (check.py v1.0.0 / gwaslab.process.py v1.4.8)
- **New tool** `gwas_process.check.py` (v1.0.0): standalone Python script that parses `*.out` / `*.err` log files produced by the staged gwaslab pipeline and prints a per-stage summary table with key QC metrics, warning/error counts, and overall pass/fail status.
- Supports checking a single study (`python gwaslab.process.check.py GWAS_ID [log_dir]`), all studies in a directory (`--all`), or only studies with problems (`--errors-only`).
- Parses: variant counts from preprocess/normalize, liftover status, chr-split count, checkref match rate + flipped/unmatched variant totals aggregated across chromosomes, and merge combined/QC-pass/COJO variant counts.
- Distinguishes real errors (Traceback, `*Error:`, `Illegal instruction`, `[ERROR]`) from known-benign upstream warnings (gwaslab FutureWarning/UserWarning/SettingWithCopyWarning, matplotlib, htslib `[W::]`).
- 📝**Updated**: `README.md` — added `gwas_process.check.py` to the HPC files table and added a dedicated `🩺 Check run status` section with usage examples and sample output.
- 🛠️**Updated**: Bumped `gwas_process.py` to v1.4.8 (`2026-03-18`) to mark this as a versioned release.

## 2026-03-18 🔇 Suppress FutureWarning in run_merge pd.concat (v1.4.7)
- **Warning**: `FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated` was emitted from line 903 during the merge stage. Root cause: per-chromosome shards for continuous traits (no `N_cases`/`N_controls`) contain all-NA entries in those nullable integer columns. When `pd.concat` sees a mix of all-NA and populated shards it warns about dtype inference, even though the parquet-preserved `Int64` dtypes are already correct and consistent across all shards.
- 🐛**Fixed** `gwas_process.py` — `run_merge()` now (1) filters out genuinely empty DataFrames before concat as a safety guard, and (2) wraps `pd.concat` in a `warnings.catch_warnings()` context that suppresses only this specific FutureWarning. The current concat behaviour is exactly correct for our use case; the suppression will be revisited if pandas changes its dtype inference in a way that affects results.
- 🛠️**Updated**: Added `import warnings` to the import block. Bumped version to `1.4.7` (`2026-03-18`).

## 2026-03-18 🐛 SIGILL crash on older HPC compute nodes — polars requires AVX2 (v1.4.6)
- **Bug**: All per-chromosome SLURM array jobs crashed with `Illegal instruction (core dumped)` preceded by a `polars` RuntimeWarning: `Missing required CPU features: avx2, fma, bmi1, bmi2, lzcnt, movbe`. The standard `polars` wheel on PyPI is compiled with AVX2/FMA intrinsics; older HPC compute nodes (pre-Haswell microarchitecture) lack those instruction-set extensions. Setting `POLARS_SKIP_CPU_CHECK=1` only suppresses the Python-level RuntimeWarning — the binary still faults the moment any AVX2 instruction executes. The crash occurs inside `gwaslab` which imports `polars` internally.
- 🐛**Fixed** `environment.yml` — replaced `"polars>=1.27.0"` with `"polars-lts-cpu>=1.27.0"`. The `polars-lts-cpu` PyPI package is the official CPU-compatible build of Polars, compiled for SSE2/SSE4 without AVX2 or FMA requirements. It provides an identical public API and can run on any x86-64 node regardless of CPU generation. The tradeoff is a modest performance reduction on modern nodes (typically 10–20 % slower for vectorised operations), which is negligible compared to the I/O and Python overhead in this pipeline.
- ⚠️ **Action required**: rebuild the conda environment after pulling this change: `mamba env remove -n gwas2cojo && mamba env create -f environment.yml`. If the environment must be updated in-place without rebuilding: `pip uninstall polars && pip install "polars-lts-cpu>=1.27.0"`.
- 🛠️**Updated**: Bumped version to `1.4.6` (`2026-03-18`).

## 2026-03-18 🐛 Categorical re-encoding by gwaslab __init__ breaks flip_allele_stats in per-chr stages (v1.4.5)
- **Bug**: `process-check-ref` (and downstream per-chromosome stages) crashed with `TypeError: Cannot setitem on a Categorical with a new category, set the categories first` inside `gwaslab.flip_allele_stats()`, even though `load_chrom_parquet()` already converted Categorical columns to `object` before passing the DataFrame to `make_sumstats_from_chrom_df()`. Root cause: gwaslab's `gl.Sumstats.__init__()` calls `basic_check()` internally, which re-encodes `EA` and `NEA` as `pd.Categorical`. Because the input is a per-chromosome shard, each column's category set only contains the allele values observed on that chromosome. When `flip_allele_stats()` then tries to swap `EA↔NEA` for 117,715 variants (e.g. LAS chromosome 9), it attempts to assign an `EA` value that is present in `NEA`'s category set but absent from `EA`'s subset, causing the pandas error. The earlier `v1.4.2` fix in `load_chrom_parquet()` was not sufficient because it converted before construction, and construction undoes it.
- 🐛**Fixed** `gwas_process.py` — `make_sumstats_from_chrom_df()` now runs a second Categorical-to-object conversion on `gwas_obj.data` immediately *after* `make_sumstats_object()` returns. All columns identified by `select_dtypes(include="category")` (typically `EA`, `NEA`, `SNPID`) are converted to plain `object` dtype. This conversion is applied once at object-creation time and persists through all downstream per-chromosome processing steps (`check_ref`, `infer_strand`, `assign_rsid`, `check_af`). The `load_chrom_parquet()` conversion is kept as a pre-construction safety net.
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
- 🐛**Fixed** `gwas_process.array_for_submit.sh` — the conf-loading stanza now checks the environment variable `GWAS2COJO_CONF` first (exported by the submit scripts, which run on the login node and always have the correct absolute path). The `BASH_SOURCE`-relative lookup is retained as a fallback for direct local invocation only. An improved error message names all three possible causes when the conf is still not found.
- 🛠️**Updated** `gwas_process.submit_staged.sh`, `gwas_process.submit.sh` — both scripts now `export GWAS2COJO_CONF="${CONF}"` immediately after sourcing the conf. SLURM propagates all exported environment variables to job environments by default (`--export=ALL`), so the absolute path is reliably available inside every job regardless of which spool directory SLURM uses.

## 2026-03-18 🔒 Removed hardcoded site-specific paths; added gwas2cojo.conf
- 🛠️**Updated**: Bumped version to `1.4.3` (`2026-03-18`).
- 🔒**Removed** all hardcoded HPC-specific paths and institutional email addresses (`@umcutrecht.nl`) from every tracked file in the repository so the codebase is clean for public use.
- 🆕**Added** `gwas2cojo.conf.example` — a single site-configuration template containing five variables (`PYTHON_SCRIPT`, `REF_DIR`, `OUT_BASE`, `CONDA_ENV`, `EMAIL`). Users copy it to `gwas2cojo.conf` (gitignored) and fill in their local values once.
- 🛠️**Updated** `gwas_process.array_for_submit.sh`, `gwas_process.submit.sh`, `gwas_process.submit_staged.sh`, `gwas_process.cleanup.sh` — replaced per-script `USER CONFIGURATION` blocks with a uniform conf-loading stanza (`source "${SCRIPT_DIR}/gwas2cojo.conf"`). All four scripts now emit a clear error with instructions if `gwas2cojo.conf` is missing.
- 🛠️**Updated** `.gitignore` — added `gwas2cojo.conf` and `gwas_list.txt` so local site settings and study lists are never accidentally committed.
- 🛠️**Updated** `gwaslab.download_refs.py` — replaced hardcoded `DEFAULT_REF_DIR` path and docstring with placeholder values.
- 🆕**Added** `gwas_list.example.txt` — a minimal three-study example config (`CAD_Aragam`, `CHARGE_CAC_EA`, `AF`) with placeholder `/path/to/gwas_datasets/` prefixes, HEADER comments, and resource annotations. Serves as the committed template; users copy to `gwas_list.txt` (gitignored) and update paths.
- 🛠️**Updated** `gwas2cojo.py`, `gwas2cojo-verify.py` — replaced `@umcutrecht.nl` banner addresses with obfuscated personal addresses (`lennart[at]landsmeer[dot]email`, `s.w.vanderlaan[at]gmail[dot]com`), matching the format already used in `gwas_process.py` and the README licence block.
- 🛠️**Updated** `README.md` — added a `⚙️ One-time site setup` section explaining the `gwas2cojo.conf` workflow; added `gwas2cojo.conf.example` and `gwas_list.example.txt` to the HPC files table; replaced remaining HPC paths and institutional emails throughout.

## 2026-03-18 🐛 Categorical dtype crash in per-chr stages and ZeroDivisionError in plot_daf
- 🛠️**Updated**: Bumped version to `1.4.2` (`2026-03-18`).
- **Bug 1**: All per-chromosome `process-check-ref` and `process-infer-strand` jobs failed with `TypeError: Cannot setitem on a Categorical with a new category, set the categories first` inside gwaslab's `flip_allele_stats()`. gwaslab encodes EA/NEA/SNPID as `pd.Categorical` after `basic_check()` for memory efficiency, and parquet round-trips preserve that dtype. A per-chromosome shard's Categorical column only contains the allele categories actually present on that chromosome. When `flip_allele_stats` tries to swap EA↔NEA for a variant whose allele (e.g. an indel sequence) exists in EA's category set but not NEA's on that chromosome, pandas refuses the assignment. In the whole-genome path this never surfaced because the full Categorical across all chromosomes includes all values in both columns simultaneously.
- 🐛**Fixed** `gwas_process.py` — `load_chrom_parquet()` now converts all `pd.CategoricalDtype` columns to plain `object` dtype immediately after reading the parquet (`df.select_dtypes(include="category")`), before the DataFrame is passed to `make_sumstats_from_chrom_df()`. This only affects EA/NEA/SNPID-style string columns; numeric columns (STATUS `Int64`, CHR, POS, BETA, SE, P, N, etc.) are not Categorical and are completely unaffected. gwaslab operates identically on `object` dtype allele strings for all harmonise/check operations; Categorical is purely a memory optimisation that is not required for correctness.
- **Bug 2**: The `merge` stage failed for CHARGE_CAC_EA_AA with `ZeroDivisionError: division by zero` inside gwaslab's `plot_daf()` (`num / len(sumstats)` where `len(sumstats) == 0`). The study had too few variants with a valid DAF value after processing (EAF largely absent or all-NaN), leaving an empty subset after DAF filtering inside gwaslab's plot routine.
- 🐛**Fixed** `gwas_process.py` — wrapped `gwas_obj.plot_daf()` in both `plot_full_dataset()` and `plot_qc_dataset()` with `try/except ZeroDivisionError`. When triggered, a `logging.warning` is emitted and the DAF plot is skipped; the rest of the merge stage (Manhattan, QQ, QC, leads, COJO) continues normally.
- 🛠️**Updated** `gwas_process.submit.sh` — added `NODES`, `CPUS`, `EMAIL`, and `MAIL_TYPE` variables and passed `--nodes`, `--cpus-per-task`, `--mail-type`, `--mail-user` to the `sbatch` call. Previously these were absent, relying on the now-removed `#SBATCH` directives in `array_for_submit.sh`.
- 🛠️**Updated** `gwas_process.array_for_submit.sh` — removed `#SBATCH --mail-type=END,FAIL` and `#SBATCH --mail-user` directives from the worker script header. SLURM merges `#SBATCH` directives from the script with command-line flags rather than letting the command line override them, so the hardcoded `END` in the script was causing end-of-job emails despite both submit scripts setting `--mail-type=FAIL`. Mail settings are now solely controlled by the calling submit script. Updated the comment to correctly name both `submit.sh` and `submit_staged.sh` as the controlling scripts.

## 2026-03-18 🐛 Missing --nodes, --cpus-per-task, and --mail-* in gwaslab.process.submit_staged.sh
- 🛠️**Updated**: Bumped version to `1.4.1` (`2026-03-18`).
- **Bug**: The v1.4.0 per-chromosome refactor of `gwas_process.submit_staged.sh` dropped three SLURM job settings that were present in the earlier script: `--nodes`, `--cpus-per-task`, and `--mail-type`/`--mail-user`. As a result all submitted jobs would inherit SLURM defaults (typically 1 CPU, which starves the multi-threaded Python worker that requests `--threads 8` via `WORKER_FLAGS`), and no failure-notification emails would be sent.
- 🐛**Fixed** `gwas_process.submit_staged.sh` — added four variables to the USER CONFIGURATION block (`NODES=1`, `CPUS=8`, `EMAIL`, `MAIL_TYPE="FAIL"`) and passed `--nodes`, `--cpus-per-task`, `--mail-type`, `--mail-user` to all eight `sbatch` calls (preprocess, normalize, split, check-ref, infer-strand, assign-rsid, check-af, merge). `CPUS` is intentionally kept in sync with the `--threads N` value in `WORKER_FLAGS`.

## 2026-03-18 🐛 Wrong build passed to reference-checking stages after liftover
- **Bug**: In all staged pipeline paths (`process-check-ref`, `process-infer-strand`, `process-assign-rsid`, `process-check-af`, `qc`, and the new `merge`), `normalise_build(REFERENCE)` was used instead of `build_num` when selecting the reference FASTA, dbSNP VCF, and Sumstats build, and when setting the chromosome map for Manhattan/QQ plots. `REFERENCE` is set from `args.build` (the original input build, e.g. `"19"`) and is never updated between staged invocations. `build_num` is correctly set to `"38"` at startup for any hg19/hg18+liftover study.
- **Impact**: For any study submitted with `--liftover`, all four heavy stages and both plot-generating stages would use:
  - `hg19.fa.gz` instead of `hg38.fa.gz` in `check_ref` → coordinates are hg38, FASTA is hg19 → nearly all variants incorrectly flagged as MISREF and lost.
  - `GCF_000001405.25.gz` (hg19 dbSNP) instead of `GCF_000001405.40.gz` (hg38 dbSNP) in `assign_rsid` → rsIDs assigned from the wrong coordinate space.
  - `build="19"` in reconstructed `Sumstats` objects (per-chr and merge paths) → wrong internal build attribute for all downstream gwaslab operations.
  - `build=reference` in `plot_mqq` → hg19 chromosome-length map applied to hg38 positions → distorted Manhattan plots.
- **Why not seen before**: The staged whole-genome path (pre-v1.4.0) always OOM'd inside `process-check-ref` or later, so these stages never produced output. The per-chromosome refactor (v1.4.0) is specifically designed to make these stages complete — meaning the wrong results would be written and stored for the first time.
- 🐛**Fixed** `gwas_process.py` — replaced `normalise_build(REFERENCE)` / `REFERENCE` with `build_num` at nine call sites across `main()` and `run_merge()`:
  - `process-check-ref` per-chr and whole-genome: `run_check_ref(gwas_obj, build_num, args.ref)` (FASTA path)
  - `process-assign-rsid` per-chr and whole-genome: `run_assign_rsid(gwas_obj, build_num, args.ref, …)` (dbSNP VCF path)
  - `make_sumstats_from_chrom_df(df, build_num)` in all four per-chr stage branches and in `run_merge()`
  - `plot_full_dataset(…, build_num, …)` and `plot_qc_dataset(…, build_num, …)` in `process-check-af` (whole-genome), `qc`, and `merge`

## 2026-03-17 🆕 Per-chromosome array-job pipeline for heavy stages (v1.4.0)
- 🛠️**Updated**: Bumped version to `1.4.0` (`2026-03-17`).
- **Context**: Studies were OOM-failing at `process-check-ref`, `process-infer-strand`, `process-assign-rsid`, and `process-check-af` even at 128–256 G. The root cause is that these stages sweep large VCF files (1KG ~84 M variants, dbSNP ~1 B variants) against the full genome-wide GWAS dataset. The fix splits the dataset by chromosome before the heavy stages so each VCF-sweep job works on ~1/22 of the variants.
- 🆕**Added** `process-split` stage to `gwas_process.py` — loads `{stem}.normalize.pkl`, splits by chromosome into per-chromosome BROTLI parquets (`{stem}.chr{N}.normalize.parquet`, N = 1–26), and writes a `{stem}.chrsplit.json` manifest. CHR values follow gwaslab's `Int64` convention: 1–22 = autosomes, 23 = X, 24 = Y, 25 = nonPAR, 26 = MT. Parquets preserve the STATUS bitmask column so gwaslab state is maintained across the per-chr jobs.
- 🆕**Added** `merge` stage to `gwas_process.py` — concatenates all `{stem}.chr{N}.checkaf.parquet` shards into a single genome-wide DataFrame, recreates a gwaslab `Sumstats` object (with STATUS restored), then runs QC filtering, plots (Manhattan, QQ, DAF), lead-variant extraction, and COJO output. Replaces the separate `qc` + `cojo` stages in the per-chromosome pipeline path.
- 🆕**Added** `--chrom N` argument (int 1–26) to `gwas_process.py`. When set, `process-check-ref`, `process-infer-strand`, `process-assign-rsid`, and `process-check-af` each operate on a single chromosome shard (`{stem}.chr{N}.{prev}.parquet` → `{stem}.chr{N}.{next}.parquet`). If the shard does not exist the stage exits gracefully with exit code 0, satisfying SLURM `afterok` dependencies for the next array stage.
- 🆕**Added** helper functions: `load_chrom_parquet()`, `save_chrom_parquet()`, `make_sumstats_from_chrom_df()`, `split_by_chrom()`, `load_chrsplit_manifest()`, `run_merge()`.
- 🛠️**Updated** `gwas_process.submit_staged.sh` — the four heavy process stages are now submitted as SLURM array jobs (`--array=1-26`); a `process-split` job is inserted between `process-normalize` and the array stages; a `merge` job replaces the `qc` + `cojo` tail. Fixed resources: `process-split` 16 G / 30 min. `afterok` on an array job ID waits for all 26 tasks; absent-chromosome tasks (exit 0) satisfy the dependency automatically. Job count per study: ~107 (vs. 8 before), well within the site limit of 120,000. Updated monitor/cancel hints.
- 🛠️**Updated** `gwas_process.array_for_submit.sh` — appends `--chrom ${SLURM_ARRAY_TASK_ID}` to the Python command when running as an array task; logs the chromosome in the job header.

## 2026-03-16 🆕 Reference file download utility
- 🆕**Added**: `gwaslab.download_refs.py` — utility script and complete inventory of all gwaslab reference files, using gwaslab's built-in `gl.download_ref()` function. Active entries (AFR, EAS, AMR, SAS for both hg19 and hg38) are downloaded; all other files already present at the reference directory are listed as comments and can be uncommented to (re-)download. Covered categories: 1KG population VCFs (all six populations, hg19 + hg38), HapMap3 EAF tables, 1KG SNPID→rsID conversion tables, dbSNP v151/v157 VCFs (very large, NCBI FTP), UCSC reference FASTA, recombination maps, and Ensembl/RefSeq GTF files. The `.tbi` index is fetched automatically alongside each VCF. Target directory defaults to `/path/to/references/gwaslab/`; override with `--ref-dir`. Prints a summary via `gl.check_downloaded_ref()` on completion.
- 🛠️**Updated**: `README.md` — added a `📥 Reference file management` section with a full reference-file inventory table (keyword, filename, default status), usage instructions, and a note about Dropbox/NCBI accessibility on HPC. Added `gwaslab.download_refs.py` to the HPC helper-scripts file table.

## 2026-03-16 🛠️ Updated environment.yml and installation instructions
- 🛠️**Updated**: `environment.yml` — overhauled to reflect the full dependency set required by `gwas_process.py`. Upgraded Python from `3.11` to `3.12`. Moved all Python packages to the `pip:` block with pinned or bounded versions: `numpy>=1.21.2,<2`, `adjusttext==0.8`, `matplotlib>=3.8,<3.9`, `pandas>=1.3,!=1.5`, `pysam==0.22.1`, `scikit-allel>=1.3.5`, `scipy>=1.12`, `seaborn>=0.12`, `h5py>=3.10.0`, `pyarrow`, `polars>=1.27.0`, `sumstats-liftover==1.1.0`, `jupyter==1.0.0`, `gwaslab`, `pyliftover`, `tqdm`. `bcftools` retained as a conda dependency (bioconda channel) rather than a pip package. Replaced `defaults` channel with `nodefaults` to avoid the Anaconda commercial repository, which is not permitted at many academic institutions; all packages are sourced exclusively from `conda-forge` and `bioconda`.
- 🛠️**Updated**: `README.md` — replaced the requirements and installation sections. Now documents Python 3.12 and `bcftools` as requirements; provides two installation paths (Option A: `mamba env create -f environment.yml`; Option B: manual `mamba create` + `pip install`); updated verification command to import `gwaslab` and `polars`; updated troubleshooting guidance for dependency conflicts and bioconda `bcftools`.

## 2026-03-16 🛠️ Two-tier resource model for gwaslab.process.submit_staged.sh
- 🛠️**Updated**: `gwas_list.txt` — added two new columns: `MEM_LIGHT` (COL10) and `TIME_LIGHT` (COL11) for the moderate pipeline stages (`process-normalize`, `process-check-ref`, `qc`). The existing `MEM` (COL8) and `TIME` (COL9) columns are unchanged and continue to control the heavy VCF-sweep stages (`process-infer-strand`, `process-assign-rsid`, `process-check-af`). Note added to header: `MEM_LIGHT` should be set higher for studies with many columns or complex allele structure (e.g. the AF multi-ancestry meta-analysis required 128G at `process-check-ref` despite having fewer variants than EUR studies that passed at 64G).
- 🛠️**Updated**: `gwas_process.submit_staged.sh` — replaced per-stage fixed defaults with two script-level fallback defaults (`MEM_LIGHT_DEFAULT=64G`, `MEM_HEAVY_DEFAULT=128G`). Per-study `MEM_LIGHT`/`TIME_LIGHT` are read from COL10/COL11 and applied to all light-tier stages (`process-normalize`, `process-check-ref`, `qc`); if absent the fallbacks are used. Report table now shows both tiers alongside the job chain.
- 🛠️**Updated**: Active entries in `gwas_list.txt` — `MEM_LIGHT`/`TIME_LIGHT` assigned per study: `32G/12h` for standard EUR studies; `64G/24h` for PAN and large EUR studies; `128G/24h` for AF (known to require higher memory at `process-check-ref`).

## 2026-03-16 🆕 Fine-grained process sub-stages, staged submit, and cleanup (v1.3.0)
- 🛠️**Updated**: Bumped version to `1.3.0` (`2026-03-16`).
- 🆕**Added**: Five `--stage process-*` sub-stages to `gwas_process.py`, splitting the monolithic process stage by memory profile. Each sub-stage saves a pickle checkpoint so subsequent stages can be submitted as independent SLURM jobs with their own resources:
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
- 🆕**Added**: `gwas_process.submit_staged.sh` — new script that submits one SLURM job per stage per study with `--dependency=afterok` chaining. If a stage fails, SLURM cancels all downstream stages for that study automatically; other studies are unaffected. `MEM` and `TIME` from `gwas_list.txt` are applied to the two heaviest stages (`process-infer-strand` and `process-assign-rsid`); all other stages use fixed resource defaults defined at the top of the script. `process-assign-rsid` is omitted when `--dbsnp` is absent from `WORKER_FLAGS`.
- 🆕**Added**: `gwas_process.cleanup.sh` — removes all intermediate checkpoint files after a successful run. Final outputs (`.parquet`, `.tsv.gz`, `.qc.*`, `.cojo.gz`, `.leads.tsv`, `.log`, `PLOTS/`) are never touched. Supports `--study NAME`, `--all`, or `--config gwas_list.txt` scope; `--dry-run` prints what would be deleted without removing; `--keep-raw-pkl` preserves the final raw pickle; `--remove-qc-pkl` also removes the QC pickle (kept by default).

## 2026-03-16 🆕 Pipeline staging in gwaslab.process.py (v1.2.0)
- 🛠️**Updated**: Bumped version to `1.2.0` (`2026-03-16`).
- 🆕**Added**: `--stage` flag to `gwas_process.py` with four stages: `preprocess`, `process`, `qc`, and `cojo` (plus `all`, the default, which preserves the existing end-to-end behaviour). Each stage can be submitted as a separate SLURM job with its own `--mem` and `--time`, allowing resource-light stages to run at `64G / 48h` while memory-intensive steps (`process`: `check_ref`, `infer_strand2`, `assign_rsid`) can be given `128G–256G / 96h` independently.
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
- 🆕**Added**: `--no-pickle` flag to `gwas_process.py`. When set, `.pkl` files are skipped for both raw and QC outputs, reducing peak memory and disk usage on the save step. The gwaslab `.log` file is still written regardless. Note: `--only-qc` requires a pickle from a prior run, so it is incompatible with `--no-pickle`.
- 🛠️**Updated**: `save_raw_outputs()` and `save_qc_outputs()` — eliminated the parquet read-back pattern (`pd.read_parquet(parquet_path)`) that was used to generate the TSV.GZ. Both functions now write the TSV directly from the in-memory `gwas_obj.data` / `gwas_obj_qc.data`, avoiding a full extra copy of the data just to write one file.
- 🛠️**Updated**: `main()` — added `del gwas_data; gc.collect()` immediately after `plot_raw_histograms()` (the last use of the raw DataFrame). This frees the raw pandas DataFrame before the heavy processing steps (`check_ref`, `infer_strand2`, `assign_rsid`, `check_af2`), preventing two full-size DataFrames from coexisting in RAM throughout the pipeline.
- 🛠️**Updated**: `main()` — added `del gwas_obj; gc.collect()` immediately after `apply_qc()` returns `gwas_obj_qc`. The unfiltered object is freed before saving QC outputs and generating QC plots, so only one copy of the data is in memory at a time during the QC stage.
- 🛠️**Updated**: `write_cojo()` — removed unnecessary `.copy()` call (`df = gwas_obj.data.copy()` → `df = gwas_obj.data`). All downstream accesses are read-only (column selection, `astype`, constructing a new `pd.DataFrame`), so the copy was wasted memory.

## 2026-03-15 🛠️ Overhaul of SLURM submission and GWAS list
- 🛠️**Updated**: The `gwas_list.txt` file to use semicolons (`;`) as the field delimiter instead of tabs, avoiding parsing issues when paths or values contain whitespace.
- 🆕**Added**: Two new columns to `gwas_list.txt`: `MEM` (COL8, SLURM memory per job, e.g. `64G` or `128G`) and `TIME` (COL9, SLURM time limit per job, e.g. `48:00:00`), allowing resource requirements to be set individually per dataset.
- 🛠️**Updated**: `gwas_process.submit.sh` to submit one independent SLURM job per dataset instead of a single array job. Memory (`--mem`) and time (`--time`) are now read from the config file and passed to each `sbatch` call individually, so datasets with different resource needs no longer share a single limit. Each job receives its own `--job-name`, `--output`, and `--error` derived from the dataset name.
- 🛠️**Updated**: `gwas_process.array_for_submit.sh` to act as a single-dataset worker script. Removed array job logic (`SLURM_ARRAY_TASK_ID`), removed fixed `--mem`, `--time`, `--output`, and `--error` SBATCH directives (these are now set dynamically by `gwas_process.submit.sh`). The script now accepts a semicolon-delimited config line as its first argument and parses it directly.

## 2025-03-12 🛠️ Updates to GWAS list
- 🆕**Added**: New GWAS datasets to the `gwas_list.txt` file, including:
    - AFGen Roselli 2018 dataset for allele frequencies (AF) with b38 positions.
    - GLGC Graham 2021 datasets for HDL, LDL, TC, TG, and non-HDL traits in European populations.
- 🧰**Fixed**: Issue with time of the SLURM job in `gwas_process.array_for_submit.sh` to allow for longer processing times, especially for larger GWAS datasets. Updated the time limit from 1 hour to 4 hours to accommodate the increased computational demands of processing multiple large GWAS datasets.

## 2025-03-12 🛠️ Updates to GWAS list
- 🆕**Added**: New GWAS datasets to the `gwas_list.txt` file, including:
    - ISGC GigaStroke datasets for ALLSTROKE, IS, CES, LAS, and SVD subtypes.
    - CHARGE cIMT (Franceschini 2018) and CHARGE Plaque (Franceschini 2018) datasets.
- 🛠️**Updated**: The `gwas_list.txt` file to ensure consistency in formatting and correct file paths.
- 🛠️**Updated**: Changed the SLURM parameters for `gwas_process.array_for_submit.sh`. 

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
    - Updated plotting functions in `gwas_process.py` to include verbose logging and ensure that plots are saved with the correct DPI settings.
    - Handles the case where a pickle file was created and the --only-qc flag is used to regenerate plots without re-running the full pipeline.
- 🛠️**Updated**: The `LICENSE` file to correct the copyright year.
- 🛠️**Updated**: The `.gitignore` file to include new directories and files that should be ignored by git.
- 🛠️**Updated**: The `CHANGES.md` file to document the new functions and updates made to the codebase.
- 🛠️**Updated**: The `README.md` file to reflect the new functionality and provide instructions for using the new script and notebook.
- 🛠️**Updated**: The `gwas_process.py` file to include the new script for processing GWAS summary statistics and to ensure that the `stem` variable is defined in all relevant branches of the code.
- 🆕**Added**: Scripts for submitting GWAS processing jobs:
    - `gwas_process.submit.sh`: A shell script to submit a GWAS processing job to a cluster using `sbatch`.
    - `gwas_process.array_for_submit.sh`: A shell script to submit an array of GWAS processing jobs for multiple datasets or parameters. This is controlled by the `gwas_process.submit.sh` script, which can be configured to run multiple instances of the processing script with different arguments.
    - `gwas_list.txt`: A text file containing a list of GWAS datasets to be processed. This file is used by the `gwas_process.array_for_submit.sh` script to determine which datasets to process in the array job. Each line in the file should specify a GWAS dataset, and the processing script will read this file to know which datasets to run on.

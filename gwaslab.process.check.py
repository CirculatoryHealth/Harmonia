#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# gwaslab.process.check.py — Run-status checker for gwas2cojo pipeline logs
# ============================================================
# Parses *.out / *.err files written by the staged gwaslab pipeline and prints
# a per-stage summary table with key QC metrics, warning counts, and any real
# errors detected.
#
# Usage
# -----
#   # Check one study in the current directory:
#   python gwaslab.process.check.py CAD_Aragam
#
#   # Check all studies in a specific log directory:
#   python gwaslab.process.check.py --all /path/to/logs
#
#   # Only show studies with errors or missing stages:
#   python gwaslab.process.check.py --all --errors-only
#
# Log file naming convention (set by submit scripts):
#   {GWAS}_{step}_{stage}_{jobid}.{out|err}            non-array stages
#   {GWAS}_{step}_{stage}_{jobid}_{chrom}.{out|err}    per-chr array stages
#
# ============================================================
VERSION_NAME = "gwaslab_process_check"
VERSION      = "1.0.0"
VERSION_DATE = "2026-03-18"
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

import sys
import os
import re
import glob
import argparse
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Stage definitions — order matters for the summary table
# ---------------------------------------------------------------------------
STAGES = [
    "preprocess",
    "normalize",
    "split",
    "checkref",
    "inferstrand",
    "assignrsid",
    "checkaf",
    "merge",
]

ARRAY_STAGES = {"checkref", "inferstrand", "assignrsid", "checkaf"}
N_CHR_EXPECTED = 26   # 1-22 + X + Y + MT + PAR

# ---------------------------------------------------------------------------
# Warning patterns that are known-benign gwaslab/matplotlib/htslib upstream
# noise — counted separately, not treated as errors.
# ---------------------------------------------------------------------------
_KNOWN_WARN = re.compile(
    r"FutureWarning|UserWarning|SettingWithCopyWarning|"
    r"SmallSampleWarning|RuntimeWarning|\[W::"
)

# Patterns that indicate a real failure
_ERROR_RE = [
    re.compile(r"^Traceback \(most recent call last\)", re.M),
    re.compile(r"^[A-Za-z]+Error:", re.M),
    re.compile(r"Illegal instruction|core dumped"),
    re.compile(r"\[ERROR\]"),
    re.compile(r"DUE TO TIME LIMIT|CANCELLED"),
]

# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

# Anchor on the known stage names so GWAS names with underscores parse cleanly.
_FILE_RE = re.compile(
    r"^(.+?)_(\d+)_("
    + "|".join(STAGES)
    + r")_(\d+)(?:_(\d+))?\.(\w+)$"
)


def find_files(log_dir: str, gwas: str):
    """
    Return nested dict: stage → chrom_or_None → {out: path, err: path}.
    Only files whose GWAS prefix matches *gwas* exactly.
    """
    result = defaultdict(lambda: defaultdict(dict))
    for path in glob.glob(os.path.join(log_dir, f"{gwas}_*.out")) + \
                glob.glob(os.path.join(log_dir, f"{gwas}_*.err")):
        name = os.path.basename(path)
        m = _FILE_RE.match(name)
        if not m:
            continue
        file_gwas, step, stage, jobid, chrom_str, ext = m.groups()
        if file_gwas != gwas:
            continue
        chrom = int(chrom_str) if chrom_str is not None else None
        result[stage][chrom][ext] = path
    return result


def discover_studies(log_dir: str):
    """Return sorted list of all GWAS IDs that have a preprocess .out file."""
    studies = set()
    for path in glob.glob(os.path.join(log_dir, "*_preprocess_*.out")):
        m = _FILE_RE.match(os.path.basename(path))
        if m:
            studies.add(m.group(1))
    return sorted(studies)


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _read(path):
    """Read text from a file, return empty string on any error."""
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, errors="replace") as fh:
            return fh.read()
    except OSError:
        return ""


def _first(text, pattern, group=1, default=None):
    """Search *text* for *pattern*, return first match group or default."""
    m = re.search(pattern, text)
    if not m:
        return default
    val = m.group(group)
    return val.strip() if val is not None else default


def _int(text, pattern, group=1, default=0):
    """Search *text* for *pattern*, return first match group as int or default."""
    m = re.search(pattern, text)
    if not m:
        return default
    try:
        return int(m.group(group).replace(",", ""))
    except ValueError:
        return default


def _parse_err(path):
    """Return (n_warnings, n_errors, [first error snippet])."""
    text = _read(path)
    if not text:
        return 0, 0, []
    n_warn  = len(_KNOWN_WARN.findall(text))
    errors  = []
    for pat in _ERROR_RE:
        errors.extend(pat.findall(text))
    # De-duplicate and cap
    seen = set()
    unique = []
    for e in errors:
        e = e.strip()
        if e not in seen:
            seen.add(e)
            unique.append(e)
    return n_warn, len(unique), unique[:3]


def _is_done(out_text):
    """Heuristic completion check on an .out file."""
    return bool(
        re.search(r"\[SAVE\]", out_text) or
        re.search(r"\[merge\] Stage complete", out_text) or
        re.search(r"Finished", out_text)
    )


# ---------------------------------------------------------------------------
# Per-stage metric extractors
# ---------------------------------------------------------------------------

def _metrics_preprocess(text):
    """Extract key metrics from the preprocess stage .out text."""
    parts = []
    n = _first(text, r"\[SAVE\] Preprocess parquet.*?\((\S+ variants)")
    if n:
        parts.append(n)
    build = _first(text, r"\[INFO\] Build\s+:\s+(\S+)")
    if build:
        parts.append(f"build {build}")
    return "  |  ".join(parts) if parts else ""


def _metrics_normalize(text):
    """Extract key metrics from the normalize stage .out text."""
    parts = []
    n = _first(text, r"After duplicate removal:\s*([\d,]+)\s+variants remain")
    if n:
        parts.append(f"{n} after dedup")
    if "Liftover skipped" in text:
        parts.append("liftover: skipped")
    elif "Liftover complete" in text:
        ref = _first(text, r"REFERENCE updated to '(\S+)'")
        parts.append(f"liftover → {ref}" if ref else "liftover: done")
    return "  |  ".join(parts) if parts else ""


def _metrics_split(text):
    """Extract key metrics from the split stage .out text."""
    n_chr = _int(text, r"\[SAVE\] Chr-split manifest.*?\((\d+) chromosomes\)")
    if n_chr:
        return f"{n_chr} chromosomes"
    # Fallback: count individual [chrN] Saved lines
    saved = len(re.findall(r"\[chr\d+\] Saved", text))
    return f"{saved} chromosomes" if saved else ""


def _metrics_checkref_agg(chr_texts):
    """Aggregate checkref metrics across all chromosome .out texts."""
    match_rates, flipped_total, unmatched_total = [], 0, 0
    for text in chr_texts:
        r = _first(text, r"Raw Matching rate\s*:\s*([\d.]+)%")
        if r:
            try:
                match_rates.append(float(r))
            except ValueError:
                pass
        flipped_total   += _int(text, r"Variants flipped\s*:\s*([\d,]+)")
        unmatched_total += _int(text, r"Variants not on given reference sequence\s*:\s*([\d,]+)")
    if not match_rates:
        return ""
    avg = sum(match_rates) / len(match_rates)
    parts = [f"match {avg:.1f}%", f"flipped {flipped_total:,}"]
    if unmatched_total:
        parts.append(f"unmatched {unmatched_total:,}")
    return "  |  ".join(parts)


def _metrics_merge(text):
    """Extract key metrics from the merge stage .out text."""
    parts = []
    combined = _first(text, r"\[merge\] Combined:\s*([\d,]+)\s+variants")
    if combined:
        parts.append(f"combined {combined}")
    qc_n = _first(text, r"After QC[^:]*:\s*([\d,]+)\s+variants?")
    if qc_n:
        parts.append(f"QC-pass {qc_n}")
    cojo = _first(text, r"\[SAVE\] COJO\s+→.*?\(([\d,]+)\s+variants?\)")
    if cojo:
        parts.append(f"COJO {cojo}")
    return "  |  ".join(parts) if parts else ""


# ---------------------------------------------------------------------------
# Study checker
# ---------------------------------------------------------------------------

def check_study(gwas: str, log_dir: str, errors_only: bool = False):
    """Check all stages for a single GWAS study and print a summary table."""
    print(f"\n===== Checking study '{gwas}' =====")
    files = find_files(log_dir, gwas)

    if not files:
        if not errors_only:
            print(f"  No log files found for '{gwas}' in {log_dir}\n")
        return

    # -- Extract header info from preprocess .out -------------------------
    pre_text = _read((files.get("preprocess") or {}).get(None, {}).get("out", ""))
    pop   = _first(pre_text, r"Population\s+:\s+(\S+)") or "?"
    build = _first(pre_text, r"Build\s+:\s+(\S+)")       or "?"

    # -- Build per-stage rows ---------------------------------------------
    rows = []        # (stage, status_str, metric_str, n_warn, n_err, [snippets])
    any_error = False

    for stage in STAGES:
        if stage not in files:
            rows.append((stage, "— missing", "", 0, 0, []))
            any_error = True
            continue

        chroms   = sorted(files[stage].keys())
        is_array = stage in ARRAY_STAGES

        if is_array:
            n_total  = len(chroms)
            n_done   = 0
            tot_warn = 0
            tot_err  = 0
            snips    = []
            chr_texts = []

            for chrom in chroms:
                out_text  = _read(files[stage][chrom].get("out", ""))
                w, e, s   = _parse_err(files[stage][chrom].get("err", ""))
                tot_warn += w
                tot_err  += e
                snips    += s
                if _is_done(out_text):
                    n_done += 1
                chr_texts.append(out_text)

            if n_done == n_total:
                status = f"✓ {n_done}/{n_total} chr"
            elif n_done == 0:
                status = f"✗ 0/{n_total} chr"
                any_error = True
            else:
                status = f"⚠ {n_done}/{n_total} chr"
                any_error = True

            if stage == "checkref":
                metric = _metrics_checkref_agg(chr_texts)
            else:
                metric = f"{n_done}/{n_total} complete"

        else:
            chrom    = None
            out_text = _read(files[stage][None].get("out", ""))
            w, e, s  = _parse_err(files[stage][None].get("err", ""))
            tot_warn, tot_err, snips = w, e, s

            done = _is_done(out_text)
            if done:
                status = "✓ done"
            else:
                # File exists but no completion marker → still running or failed
                if files[stage][None].get("out") and os.path.isfile(
                        files[stage][None]["out"]):
                    status = "⚠ incomplete"
                    any_error = True
                else:
                    status = "— missing"
                    any_error = True

            if stage == "preprocess":
                metric = _metrics_preprocess(out_text)
            elif stage == "normalize":
                metric = _metrics_normalize(out_text)
            elif stage == "split":
                metric = _metrics_split(out_text)
            elif stage == "merge":
                metric = _metrics_merge(out_text)
            else:
                metric = ""

        if tot_err:
            any_error = True

        rows.append((stage, status, metric, tot_warn, tot_err, snips[:2]))

    # -- Print table -------------------------------------------------------
    if errors_only and not any_error:
        return

    W = 100
    print("═" * W)
    print(f"  Study: {gwas}  |  Population: {pop}  |  Build: {build}")
    print("═" * W)
    print(f"  {'Stage':<16}  {'Status':<16}  {'Key metric':<46}  {'Warn':>5}  {'Err':>4}")
    print("─" * W)

    total_warn = total_err = 0
    for stage, status, metric, nw, ne, snips in rows:
        flag = "  ◄ ERRORS" if ne else ""
        metric_str = (metric[:44] + "…") if len(metric) > 45 else metric
        print(f"  {stage:<16}  {status:<16}  {metric_str:<46}  {nw:>5}  {ne:>4}{flag}")
        for snip in snips:
            print(f"      ⚠  {snip[:88]}")
        total_warn += nw
        total_err  += ne

    print("─" * W)
    overall = "✓ COMPLETE" if not any_error else "✗ ISSUES FOUND"
    print(f"  Overall: {overall:<20}  "
          f"Warnings (gwaslab upstream): {total_warn}   Errors: {total_err}")
    print("═" * W)
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Main entry point: parse args and check studies."""
    ap = argparse.ArgumentParser(
        description="Quick run-status checker for gwaslab.process pipeline log files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gwaslab.process.check.py CAD_Aragam
  python gwaslab.process.check.py CAD_Aragam /path/to/logs
  python gwaslab.process.check.py --all
  python gwaslab.process.check.py --all /path/to/logs --errors-only
  python gwaslab.process.check.py --all | tee gwaslab.process.check_$(date +%Y%m%d_%H%M%S).log
""",
    )
    ap.add_argument("gwas",    nargs="?",
                    help="GWAS study ID to check. Ignored if --all is set.")
    ap.add_argument("log_dir", nargs="?",
                    default=".", help="Directory containing log files (default: .). Ignored if --all is set.")
    ap.add_argument("--all",   action="store_true",
                    help="Check every study found in log_dir. Overrides positional gwas argument.")
    ap.add_argument("--errors-only", action="store_true",
                    help="Only print studies that have errors or missing stages. Useful for quick scans.")
    args = ap.parse_args()

    log_dir = os.path.abspath(args.log_dir or ".")

    if args.all:
        studies = discover_studies(log_dir)
        if not studies:
            print(f"No studies found in {log_dir}")
            sys.exit(1)
        print(f"Checking {len(studies)} studies in {log_dir}\n")
        for gwas in studies:
            check_study(gwas, log_dir, errors_only=args.errors_only)
    elif args.gwas:
        check_study(args.gwas, log_dir, errors_only=args.errors_only)
    else:
        ap.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

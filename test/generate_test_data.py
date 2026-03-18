#!/usr/bin/env python3
"""
generate_test_data.py — Synthesise three GWAS summary-statistics files and a
matching genetic reference for use with gwas2cojo.py tests.

Three datasets with intentionally different column naming conventions are
generated so both the auto-detection and manual-override paths of gwas2cojo.py
can be exercised:

  Dataset 1  gwas_metal_tab.txt.gz   METAL-style  — all columns auto-detected
  Dataset 2  gwas_plink2.txt.gz      PLINK2-style — exercises A1FREQ / OBS_CT aliases
  Dataset 3  gwas_saige.txt.gz       SAIGE-style  — requires --gwas:effect/other/freq

All files are hg19, tab-delimited, gzip-compressed.  The companion run.sh
shows the exact gwas2cojo.py invocation for each dataset.

Allele flips and complement-strand variants are introduced deliberately so the
FLIP and complement-handling logic is exercised in addition to the NOP path.

Requirements: Python 3.6+, standard library only (no pandas / numpy / scipy).

Usage
-----
  python3 test/generate_test_data.py              # 10 000 variants, seed 42
  python3 test/generate_test_data.py --n 5000     # fewer variants
  python3 test/generate_test_data.py --n 10000 --seed 99
"""

import argparse
import gzip
import math
import os
import random
import sys

# ── Allele table ──────────────────────────────────────────────────────────────
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Non-palindromic (non-ambivalent) and palindromic (ambivalent) SNP pairs.
# ref = non-effect allele, alt = effect allele.
NON_AMBIV = [
    ('A', 'C'), ('A', 'G'), ('T', 'C'), ('T', 'G'),
    ('C', 'A'), ('C', 'T'), ('G', 'A'), ('G', 'T'),
]
AMBIV = [('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')]

# ── Approximate hg19 chromosome sizes (Mb) ────────────────────────────────────
CHR_LEN_MB = {
     1: 249,  2: 243,  3: 198,  4: 191,  5: 181,  6: 171,  7: 159,  8: 146,
     9: 141, 10: 136, 11: 135, 12: 134, 13: 115, 14: 107, 15: 103, 16:  90,
    17:  81, 18:  78, 19:  59, 20:  63, 21:  48, 22:  51,
}

# ── Stdlib-only random helpers ────────────────────────────────────────────────

def _normal(mu=0.0, sigma=1.0):
    """Box-Muller normal sample."""
    u1 = max(random.random(), 1e-15)
    u2 = random.random()
    return mu + sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def _beta(a=0.5, b=0.5):
    """Beta(a, b) sample via Gamma ratio."""
    x = random.gammavariate(a, 1.0)
    y = random.gammavariate(b, 1.0)
    return x / (x + y + 1e-15)


def _p_from_z(z):
    """Two-sided p-value from Z score (Abramowitz & Stegun rational approximation)."""
    z = abs(z)
    if z > 38:
        return 1e-300
    t = 1.0 / (1.0 + 0.2316419 * z)
    poly = t * (0.319381530 + t * (-0.356563782
           + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return min(1.0, 2.0 * math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi) * poly)


# ── Chromosome-weighted sampler ───────────────────────────────────────────────

def _make_chrom_sampler():
    """Return a callable that samples chromosomes weighted by length."""
    chroms = list(CHR_LEN_MB.keys())
    total  = sum(CHR_LEN_MB.values())
    cum = []
    s = 0.0
    for c in chroms:
        s += CHR_LEN_MB[c] / total
        cum.append((s, c))

    def sample():
        r = random.random()
        for threshold, c in cum:
            if r <= threshold:
                return c
        return chroms[-1]

    return sample


# ── Variant generation ────────────────────────────────────────────────────────

def generate_variants(n, seed=42):
    """
    Generate *n* synthetic biallelic SNPs with realistic allele frequencies,
    effect sizes, and p-values.

    Returns a list of dicts with keys:
        chrom, pos, rsid, ref (non-effect), alt (effect), af, beta, se, p, n_obs
    """
    random.seed(seed)
    pick_chr = _make_chrom_sampler()
    used_pos = {}   # chrom -> set(positions already chosen)
    variants = []

    for i in range(n):
        chrom = pick_chr()
        used  = used_pos.setdefault(chrom, set())

        # Pick a position not already used on this chromosome.
        for _ in range(500):
            pos = random.randint(500_000, CHR_LEN_MB[chrom] * 1_000_000 - 500_000)
            if pos not in used:
                used.add(pos)
                break

        # ~10 % palindromic (ambivalent) SNPs.
        ref, alt = random.choice(AMBIV if random.random() < 0.10 else NON_AMBIV)

        # Allele frequency drawn from Beta(0.5, 0.5) — U-shaped, realistic.
        af = round(max(0.01, min(0.99, _beta(0.5, 0.5))), 4)

        # Sample size and SE.
        n_obs = random.randint(5_000, 100_000)
        se    = max(1.0 / math.sqrt(2.0 * n_obs * af * (1.0 - af) + 1e-9), 5e-4)

        # ~5 % variants have a real (non-null) signal.
        beta = _normal(0.0, 0.12) if random.random() < 0.05 else _normal(0.0, 0.008)
        p    = _p_from_z(beta / se)

        variants.append({
            'chrom': chrom,
            'pos':   pos,
            'rsid':  f'rs{1_000_000 + i}',
            'ref':   ref,           # non-effect allele
            'alt':   alt,           # effect allele
            'af':    af,            # frequency of effect allele (alt)
            'beta':  round(beta, 6),
            'se':    round(se,   6),
            'p':     p,
            'n_obs': n_obs,
        })

    variants.sort(key=lambda v: (v['chrom'], v['pos']))
    return variants


# ── File writers ──────────────────────────────────────────────────────────────

def _gz(path):
    return gzip.open(path, 'wt', compresslevel=6)


def write_ref(variants, path):
    """
    Genetic reference for gwas2cojo.py --gen.

    Columns: chr  pos  id  a1  a2  eaf
      a1  = effect allele (ALT in VCF convention)
      a2  = non-effect allele (REF in VCF convention)
      eaf = effect-allele frequency

    gwas2cojo.py invocation requires:
      --gen:ident id  --gen:effect a1  --gen:other a2  --gen:eaf eaf
    """
    with _gz(path) as f:
        f.write('chr\tpos\tid\ta1\ta2\teaf\n')
        for v in variants:
            f.write(f"{v['chrom']}\t{v['pos']}\t{v['rsid']}\t"
                    f"{v['alt']}\t{v['ref']}\t{v['af']}\n")


def write_gwas_metal(variants, path, rng_seed=1):
    """
    Dataset 1 — METAL-style, tab-delimited.

    All columns are auto-detected by gwas2cojo.py without any --gwas:* overrides
    (except --gwas:build hg19 since the header contains no build hint).

    Columns: MarkerName  CHR  BP  Allele1  Allele2  Freq1  Effect  StdErr  P-value  N
      Allele1 = effect allele; Allele2 = non-effect allele.

    ~15 % of variants are presented on the complement strand (alleles
    complemented, frequency and beta unchanged) to exercise the translated-allele
    NOP/FLIP path in select_action().
    """
    random.seed(rng_seed)
    with _gz(path) as f:
        f.write('MarkerName\tCHR\tBP\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP-value\tN\n')
        for v in variants:
            eff, oth = v['alt'], v['ref']   # effect allele first
            if random.random() < 0.15:      # present on complement strand
                eff, oth = COMP[eff], COMP[oth]
            marker = f"{v['chrom']}:{v['pos']}:{eff}:{oth}"
            f.write(
                f"{marker}\t{v['chrom']}\t{v['pos']}\t{eff}\t{oth}\t"
                f"{v['af']:.4f}\t{v['beta']:.6f}\t{v['se']:.6f}\t"
                f"{v['p']:.3e}\t{v['n_obs']}\n"
            )


def write_gwas_plink2(variants, path, rng_seed=2):
    """
    Dataset 2 — PLINK2-style (.linear / .logistic.hybrid), tab-delimited.

    Exercises the A1FREQ and OBS_CT column aliases added in gwas2cojo.py v1.4.4.
    All columns are auto-detected without --gwas:* overrides.

    Columns: CHR  POS  ID  A1  A2  A1FREQ  BETA  SE  P  OBS_CT
      A1     = tested / effect allele
      A1FREQ = effect-allele frequency
      OBS_CT = observation count (maps to N)

    ~20 % of variants are switched (A1 = non-effect allele) with frequency and
    beta negated, exercising the allele-flip (ACT_FLIP) branch.
    """
    random.seed(rng_seed)
    with _gz(path) as f:
        f.write('CHR\tPOS\tID\tA1\tA2\tA1FREQ\tBETA\tSE\tP\tOBS_CT\n')
        for v in variants:
            a1, a2   = v['alt'], v['ref']
            freq, beta = v['af'], v['beta']
            if random.random() < 0.20:      # switch effect/non-effect allele
                a1, a2     = v['ref'], v['alt']
                freq, beta = round(1.0 - v['af'], 4), round(-v['beta'], 6)
            f.write(
                f"{v['chrom']}\t{v['pos']}\t{v['rsid']}\t{a1}\t{a2}\t"
                f"{freq:.4f}\t{beta:.6f}\t{v['se']:.6f}\t{v['p']:.3e}\t{v['n_obs']}\n"
            )


def write_gwas_saige(variants, path, rng_seed=3):
    """
    Dataset 3 — SAIGE-style, tab-delimited.

    SAIGE convention: Allele1 = non-effect (REF), Allele2 = effect (ALT).
    Because Allele1/Allele2 map to non-effect/effect (opposite of the METAL
    convention that gwas2cojo.py assumes by default), manual column overrides
    are required:
        --gwas:effect Allele2  --gwas:other Allele1  --gwas:freq AF_Allele2

    Columns: CHR  POS  SNPID  Allele1  Allele2  AF_Allele2  BETA  SE  p.value  N
      Allele1    = non-effect allele (REF)
      Allele2    = effect allele (ALT)
      AF_Allele2 = effect-allele frequency

    ~10 % of variants are presented on the complement strand.
    """
    random.seed(rng_seed)
    with _gz(path) as f:
        f.write('CHR\tPOS\tSNPID\tAllele1\tAllele2\tAF_Allele2\tBETA\tSE\tp.value\tN\n')
        for v in variants:
            a1 = v['ref']   # Allele1 = non-effect
            a2 = v['alt']   # Allele2 = effect
            if random.random() < 0.10:      # complement strand
                a1, a2 = COMP[a1], COMP[a2]
            snpid = f"{v['chrom']}:{v['pos']}:{a1}:{a2}"
            f.write(
                f"{v['chrom']}\t{v['pos']}\t{snpid}\t{a1}\t{a2}\t"
                f"{v['af']:.4f}\t{v['beta']:.6f}\t{v['se']:.6f}\t"
                f"{v['p']:.3e}\t{v['n_obs']}\n"
            )


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--n', dest='n_variants', type=int, default=10_000,
        metavar='N',
        help='Number of variants to generate (default: 10 000)',
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Master random seed (default: 42)',
    )
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(__file__))

    print(f"Generating {args.n_variants:,} synthetic variants  (seed={args.seed}) …")
    variants = generate_variants(args.n_variants, seed=args.seed)
    print(f"  Chromosomes represented : {len({v['chrom'] for v in variants})}")
    print(f"  Variants with p < 5e-8  : {sum(1 for v in variants if v['p'] < 5e-8)}")

    writers = [
        ('ref.txt.gz',            write_ref,        {}),
        ('gwas_metal_tab.txt.gz', write_gwas_metal, {'rng_seed': 1}),
        ('gwas_plink2.txt.gz',    write_gwas_plink2,{'rng_seed': 2}),
        ('gwas_saige.txt.gz',     write_gwas_saige, {'rng_seed': 3}),
    ]

    print()
    for fname, fn, kwargs in writers:
        path = os.path.join(out_dir, fname)
        fn(variants, path, **kwargs)
        size_kb = os.path.getsize(path) / 1024
        print(f"  {fname:<30}  {size_kb:6.1f} KB")

    print(f"\nDone.  Regenerate at any time:  python3 test/generate_test_data.py")
    print(f"Run tests:                      bash test/run.sh")


if __name__ == '__main__':
    main()

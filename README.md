[gwas2cojo](https://github.com/CirculatoryHealth/gwas2cojo)<img align="right" height="200" src=logo/fulllogo_transparent.png>
============
[![Languages](https://skillicons.dev/icons?i=bash,py)](https://skillicons.dev) 

`gwas2cojo.py` is a public python script that aligns a public GWAS dataset to a genetic reference, to enable large scale cross dataset comparisons with public available GWAS datasets. Among others, it tries to deal with different dataformats and different genome builds.
Most importantly, it aligns the variant notation such that swapped, translated, wrong or ambiguous ambivalent allels are corrected or removed.

It generates a [COJO] file:

```
SNP       A1  A2  freq    b       se      p       n
rs1001    A   G   0.8493  0.0024  0.0055  0.6653  129850
rs1002    C   G   0.03606 0.0034  0.0115  0.7659  129799
rs1003    A   C   0.5128  0.045   0.038   0.2319  129830
```

[COJO]: https://yanglab.westlake.edu.cn/software/smr/#Overview

# Usage

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

# References

You will need a reference to map the data to. You can create your own, or use the one we created [one based](https://blog.llandsmeer.com/1kGp3.ref.1maf.nonbia.sumstats.gz) on the 1000G phase 3 data for Europeans. This is filtered based on MAF>1% and excludes non-bi-allelic and duplicate variants.

# Requirements

The script expects Python 3.6+. 

# 🧩 Installation

We recommend using `mamba` (a faster drop-in replacement for `conda`).

## 1️⃣ Create and activate the environment

```
mamba env create -f environment.yml
mamba activate gwas2cojo
```

## 2️⃣ Verify installation

Check that the required Python modules are available:

```
python -c "import numpy, pyliftover; print('OK')"
```

If you see OK, the environment is ready.

## 3️⃣ (Optional) Test the script

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

# 🔧 Troubleshooting

## When `pyliftover` is missing

If `pyliftover` is missing, install it manually:

```
mamba install -c conda-forge pyliftover
```

## Using `conda` instead of `mamba`

If you prefer `conda`:

```
conda env create -f environment.yml
conda activate gwas2cojo
```

# License

```
The MIT License (MIT)
Copyright (c) 1979-2025 Lennart P.L. Landsmeer (lennart[at]landsmeer[dot]email) & Sander W. van der Laan (s.w.vanderlaan[at]gmail[dot]com.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the \'Software\'), to deal in the Software without restriction,       ')
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED \'AS IS\', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT   ')
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

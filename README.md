# gwas2cojo

`gwas2cojo.py` is a public python script that aligns a public GWAS dataset to a genetic reference, to enable large scale cross dataset comparisons with public available GWAS datasets. Among others, it tries to deal with different dataformats and different genome builds.
Most importantly, it aligns the variant notation such that swapped, translated, wrong or ambiguous ambivalent allels are corrected or removed.

It generates a [COJO] file:

```
SNP       A1  A2  freq    b       se      p       n
rs1001    A   G   0.8493  0.0024  0.0055  0.6653  129850
rs1002    C   G   0.03606 0.0034  0.0115  0.7659  129799
rs1003    A   C   0.5128  0.045   0.038   0.2319  129830
```

[COJO]: https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis

You'll need a reference to map the data to, we created [one based](https://blog.llandsmeer.com/1kGp3.ref.1maf.nonbia.sumstats.gz) on the 1000G phase 3 data for Europeans. This is filtered based on MAF>1% and excludes non-bi-allelic and duplicate variants.

Usage:

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

Full list of options:

```
$ python3 ./gwas2cojo.py -h
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                    QTLTools CONVERT GWAS FOR SMR


* Written by         : Lennart Landsmeer | lennart[at]landsmeer[dot]email
* Suggested for by   : Sander W. van der Laan | s.w.vanderlaan[at]gmail[dot]com
* Last update        : 2020-04-01
* Name               : gwas2cojo
* Version            : v1.4.0

* Description        : To assess pleiotropic effects using Summarized-data Mendelian Randomization (SMR)
                       of molecular QTLs on (selected) traits, summary statistics from genome-wide
                       association studies (GWAS) are converted to the GWAS-COJO format.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
usage: gwas2cojo.py [-h] [-o cojo] [-r txt] [-rr] [-g file.stats.gz] --gwas file.txt.gz. [--header-only]
                    [--fmid MID] [--fclose CLOSE] [--ignore-indels] [--gwas:effect COLUMN]
                    [--gwas:other COLUMN] [--gwas:freq COLUMN] [--gwas:beta COLUMN] [--gwas:or COLUMN]
                    [--gwas:se COLUMN] [--gwas:p COLUMN] [--gwas:chr-bp COLUMN] [--gwas:chr COLUMN]
                    [--gwas:bp COLUMN] [--gwas:build BUILDID] [--gwas:n COLUMNS] [--gwas:sep DELIMITER]
                    [--gwas:header:remove STRING] [--gwas:default:p VALUE] [--gwas:default:beta VALUE]
                    [--gwas:default:se VALUE] [--gwas:default:chr VALUE] [--gwas:default:n VALUE]
                    [--gen:ident COLUMN] [--gen:chr COLUMN] [--gen:bp COLUMN] [--gen:effect COLUMN]
                    [--gen:other COLUMN] [--gen:eaf COLUMN] [--gen:oaf COLUMN] [--gen:maf COLUMN]
                    [--gen:minor COLUMN] [--gen:build COLUMN]
gwas2cojo.py: error: the following arguments are required: --gwas
```

See also [this link] for more background and an additional explanation.

[this link]: https://blog.llandsmeer.com/tech/2019/12/28/gwas2cojo.html

# License

```
The MIT License (MIT)
Copyright (c) 1979-2023 Lennart P.L. Landsmeer & Sander W. van der Laan

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

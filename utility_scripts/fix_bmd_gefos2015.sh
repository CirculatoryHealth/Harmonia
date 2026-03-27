#!/bin/bash

echo "Processing GEFOS 2015 datasets to fix sample size (N) in column 18..."

for study in "fn2stu:49988" "fa2stu:10805" "ls2stu:44731"; do
    name="${study%%:*}"
    n="${study##*:}"
    zcat _GEFOS/2015/${name}.MAF0.005.pos.out.gz \
    | awk -F'\t' -v OFS='\t' -v n="$n" 'NR>1{$18=n} 1' \
    | gzip > _GEFOS/2015/${name}.MAF0.005.pos.fixed.out.gz
    echo "Done: ${name} → N=${n}"
done

echo ""
echo "Verifying fixed sample sizes in column 18:"

for f in fn2stu fa2stu ls2stu; do
    echo -n "${f}: "
    zcat _GEFOS/2015/${f}.MAF0.005.pos.fixed.out.gz | awk -F'\t' 'NR==2{print $18}'
done

echo ""
echo "All done! Updated files are suffixed with .fixed.out.gz"

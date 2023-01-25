set -e

python3 ../gwas2cojo.py \
    --gen A \
    --gwas B \
    --gen:ident ident \
    --gen:effect effect \
    --gen:other other \
    --gen:eaf eaf \
    -o /tmp/gwas2cojo_test \
    --output-pos | grep -v '^+'

echo
echo "# Genetic reference (--gen):"
cat A | column -t
echo
echo "# To be aligned (--gwas):"
cat B | column -t
echo
echo "# Output (-o)"
cat /tmp/gwas2cojo_test | column -t

# combine variant and annotation info into one convenient table
# magic numbers (Fst..) will be joined in R
python dump_vars.py data/lp2.sorted.gff3 data/lp2-var-filtered.vcf.gz > data/variant-table.tsv


# combine variant and annotation info into one convenient table
# magic numbers (Fst..) will be joined in R
python dump_vars.py data/lp2.sorted.gff3 data/lp2-var-filtered.vcf.gz > data/variant-table.tsv

python dump_vars.py data/lp2.sorted.gff3.gz data/lp2-var-filtered.vcf.gz | pv -l > data/variant-table.tsv

# find the longest gene models in zebra finch ensemble annotation
<taeGut1/annot/ensGenes.sorted.bed.gz zcat|mawk '{print $3 - $2;}'|sort -rn|head -20
# ensGenes max is <1M, 6 models is >500k
# refSeq has one 1036509, the second is 495217

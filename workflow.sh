# combine variant and annotation info into one convenient table
# magic numbers (Fst..) will be joined in R
python dump_vars.py data/lp2.sorted.gff3 data/lp2-var-filtered.vcf.gz > data/variant-table.tsv

python dump_vars.py data/lp2.sorted.gff3.gz data/lp2-var-filtered.vcf.gz | pv -l > data/variant-table.tsv

# find the longest gene models in zebra finch ensemble annotation
<taeGut1/annot/ensGenes.sorted.bed.gz zcat|mawk '{print $3 - $2;}'|sort -rn|head -20
# ensGenes max is <1M, 6 models is >500k
# refSeq has one 1036509, the second is 495217

ANNOTS=~/brno3/genomes/taeGut1/annot
ANN="$ANNOTS/ensGenes.bed.gz $ANNOTS/refSeqGenes.bed.gz"
zcat $ANN | 
  bedtools sort | 
  bedtools intersect -wa -sorted -a - -b <( <80-islands/islands.bed tr -d "\r" | bedtools sort ) \
> 80-islands/genes-in-islands.bed

<80-islands/genes-in-islands.bed cut -f4 | grep ENS > 80-islands/genes-in-islands.ensg
<80-islands/genes-in-islands.bed cut -f4 | grep -v ENS > 80-islands/genes-in-islands.refseq


# pick song related genes
# the paper states that blue, orange and dark green modules are 'song modules'
<data-song/neuron_10977_mmc3.txt awk '$6 == "blue" || $6 = "orange" || $6 == "darkgreen"' > data-song/song-modules.tsv
# 689 probes

# switch to MetaCentrum

# go search it in the genome, because many of the genes were not annotated back then..
module add blat-suite-34

# extract the probe sequences
<81-song/song-modules.tsv awk '{print ">seq" NR; print $3;}' > 81-song/song-modules.fa
GENOMES=~/brno3/genomes
GENOME=$GENOMES/taeGut1/taeGut1.2bit

# search! 
blat -fastMap -dots=10 $GENOME 81-song/song-modules.fa 81-song/probe-blat.psl

# check the number of hits for fastMap
<81-song/probe-blat.psl tail -n +6 | cut -f10|sort|uniq|wc -l
# 476 - not good for 689 probes

# search again, consider -fine and -noHead
blat -dots=10 $GENOME 81-song/song-modules.fa 81-song/probe-blat-full.psl
# 634 matches looks better

# dxy calculation
python dxy.py data/lp2-var-filtered.vcf.gz populations > data/dxy.tsv

# statitstics for the paper
<data-genome/lp2-var-filtered.vcf.gz zcat| grep -v '^#' | grep -v INDEL | awk '$6 > 10' | wc -l
# 300433

# melt the variants to long tabular format
<data-genome/lp2-var-filtered.vcf.gz zcat | awk '$6 > 10 || /^#/' | python melt_vcf.py > data-genome/vars-melted.tsv

#
# some stats for results section in the paper
#
samtools depth 63-smalt-all/alldup.bam > 63-smalt-all/alldup.depth
<63-smalt-all/alldup.depth mawk '{s+=$3} END{ print s }'
# 1356041442

<51-liftover-all/lp2.fasta.fai mawk '{s+=$2} END{ print s }'
# 94120451

echo "scale=2; 1356041442 / 94120451" | bc
# average coverage: 14.40 

# extract depth for all SNPs over quality 10
 <71-var-all/lp2-var-filtered.vcf.gz zcat | 
   mawk '$1 !~ /^#/ && $0 !~ /INDEL/ && $6 > 10' | 
   egrep -o 'DP=[^;]+;' | 
   sed -e 's/DP=//' -e 's/;//' \
> 71-var-all/lp2-var-filtered.depths

# count sequences assigned to zf genome
<51-liftover-all/lp2.gff3 grep -c mRNA
# 61753

# count unique genes which were matched in zf genome
<51-liftover-all/lp2.gff3 grep coords | egrep -o 'Name=[^;]*' | uniq | sort -u | wc -l
# 9332

# conunt sequences with protein protein reciprocals
<51-liftover-all/lp2.gff3 grep mvz-annot | grep -c CDS
# 6611


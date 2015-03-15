# running mvz annotate on metacentrum

# against zebra finch proteins
module add blast-2.2.26
module add blast+-2.2.27
module add bioperl
export PATH=/storage/brno2/home/liborm/sw/cd-hit-v4.6.1-2012-08-27:/storage/brno2/home/liborm/sw/exonerate-2.2.0/src/program:$PATH
export FRAMEDP=/storage/brno2/home/liborm/sw/framedp-1.2.0

# metacentrum
CPUS=$PBS_NUM_PPN

CPUS=28


cd 30-mvz-annot-tg

IN=../10-newbler
awk '/^>/ {print $1;} /^[^>]/' $IN/454Isotigs.fna > $( basename $IN ).fasta

# sample
awk '/isotig01000/ {exit;} {print;}' 10-newbler.fastax > 10-newbler-1k.fasta

# against zebra finch genome

perl /storage/brno2/home/liborm/sw/mvz-trans/6-annotationTranscriptome_without_bioperl.pl -a . -b ../../genomes/taeniopygia_guttata/pep/Taeniopygia_guttata.taeGut3.2.4.73.pep.all.fa -d /storage/brno2/home/liborm/sw/framedp-1.2.0 -e $CPUS -f tg_proteins.txt 2> stderr > stdout


# against flycatcher
# 38 minutes wall time, 3:40 cpu time on 28 CPUs, 200 MB memory
perl /storage/brno2/home/liborm/sw/mvz-trans/6-annotationTranscriptome_without_bioperl.pl -a . -b ../../genomes/ficedula_albicollis/pep/Ficedula_albicollis.FicAlb_1.4.73.pep.all.fa -d /storage/brno2/home/liborm/sw/framedp-1.2.0 -e $CPUS -f fa_proteins.txt 2> stderr > stdout


# combined = get the best hit.. (append to logs to continue broken jobs)
perl /storage/brno2/home/liborm/sw/mvz-trans/6-annotationTranscriptome_without_bioperl.pl -a $( pwd ) -b proteins-tgfa.fa -d $FRAMEDP -e $CPUS -f annots-tgfa.txt 2>> stderr >> stdout

# ~60 K transcripts takes 50 CPU hours, 1 GB memory
# comments to mvz pipeline
# alphanumeric sort for contig output - or 0 prefixed names

# count the really annotated transcripts
OUT=22-newbler/22-newbler.fasta.annotated
awk '(/^>/ && NF > 1)' $OUT | wc -l


# get rid of the sequences that were not annotated by any reference protein
seq_filter_by_id.py 22-newbler.gff3 1 22-newbler/22-newbler.fasta.annotated fasta 22-newbler.fasta.filt -

# add the mvz annotations to scaffold
ANNOTS="50-liftover-a/22-newbler.gmap_s.gff3.gz-lo.gff3 35-mvz-allsamp-tgfa/22-newbler.gff3"
# .. scaffold.py $INFILE $ANNOTS $OUT/$GNAME.fasta $OUTGFF


# realign the reads with srma
# srma is not suited for 454 .. found in the docs;)
# <chromosome name>:<start position>-<end position>

SRMA=~/brno2/sw/srma/c-code/srma
REF=51-liftover-all/lp2.fasta
awk '{print $1 ":1-" $2 + 1;}' $REF.fai | parallel -j $CPUS $SRMA -i $OUT/alldup.bam -r $REF -o $OUT/srma-{#}.bam -R {}

###
# freebayes
REF=51-liftover-all/lp2.fasta
GNAME=$( echo ${REF##*/} | cut -d. -f1 )
BAM=63-smalt-all/alldup.bam
OUT=78-freebayes
awk '{print $1;}' $REF.fai | parallel -j $CPUS "freebayes -f $REF -r {} $BAM > $OUT/${GNAME}-{}.vcf"

# join the results
OFILE=$OUT/variants-raw.vcf

# headers
FILE=$( find $OUT -name ${GNAME}-*.vcf | sort | head -1 )
< $FILE egrep '^#' > $OFILE
# the rest in .fai order
awk '{print "'$OUT/$GNAME'-" $1 ".vcf";}' $REF.fai | parallel -j 1 "egrep -v '^#' {} >> $OFILE"
# the rest, ignore order
cat $OUT/${GNAME}-*.vcf | egrep -v '^#' >> $OFILE
< $OFILE $VCFLIB/vcffilter -f 'QUAL > 20' > $OUT/variants-qual.vcf


###
## get one bam file per individual
cat 63-smalt-all/readgroups.txt | python rg2merge.py > merge.sh
. merge.sh
 OUT=63-smalt-all bash merge.sh
# one job for IO bound task (?)
parallel -j 1 samtools index {} ::: 63-smalt-all/l*.bam
# check the correct count of output files;)
awk '{print $3;}' 63-smalt-all/readgroups.txt | sort -u | wc -l

# try pyrohmm
# one bam per individual
PYRO=~/brno2/sw/pyrohmmvar/pyrohmmvar
REF=51-liftover-all/lp2.fasta
PYROCONF=~/brno2/sw/pyrohmmvar/hmm_parameter_config
OUT=72-pyro-all
parallel -j $CPUS "$PYRO -b {} -f $REF -m $PYROCONF > $OUT/{/.}.txt" ::: 63-smalt-all/l*.bam

# less stringent params
parallel -j $CPUS "$PYRO -b {} -f $REF -m $PYROCONF -M 1 -B 15 -t 30 > $OUT/{/.}.txt" ::: 63-smalt-all/l*.bam

# TODO convert pyrohmm results to vcf, merge individuals
Column 1	The chromosome name
Column 2	The genomic coordinate started from '1'
Column 3	The variant category, homozygous or heterozygous
Column 4	The reference allele
Column 5	The called allele
Column 6	The genotype posterior probability
Column 7	The variant quality score
Column 8	The median mapping quality within a window
Column 9	The median base quality within a window

# tableau importable format


#
# altered vcf for mvz filtering pipeline
#
. ~/brno2/sw/perl-local/activate
module add parallel
module add samtools-0.1.19

SCAFFOLD=51-liftover-all/lp2.fasta
ALIGNS=63-smalt-all/alldup.bam

# outputs
OUT=73-var-mvz
VARIANTS=$OUT/lp2-var
mkdir -p $OUT

vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j $CPUS "samtools mpileup -BIDSu -L 10000 -f $SCAFFOLD -r {} $ALIGNS | bcftools view -bcg - > $OUT/part-{}.bcf" | tee $OUT/spectra.txt

# the correct version, according to -outdated- .sh script from latest ngsClean git;)
OUT=77-var-mvz2
VARIANTS=$OUT/lp2-var
mkdir -p $OUT

vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j $CPUS "samtools mpileup -u -AIDS -q 0 -Q 20 -C 50 -f $SCAFFOLD -r {} $ALIGNS | bcftools view -bgI - > $OUT/part-{}.bcf" | tee $OUT/spectra.txt

# -> merge the dataset

NGSCLEAN=~/brno2/sw/ngsClean
N_IND=15
ANC_SEQ=$SCAFFOLD

# SNPcleaner.pl
# -v	process nonvariant sites (in addition to varients)
# -a	INT	minimum number of alternate alleles per site
# -k	INT	minimum number of individuals with at least [-u INT]X coverage (requires SNPcleaner -u and mpileup -D)
# -u	INT	minimum individual coverage threshold used for -k (requires SNPcleaner -k and mpileup -D) 
# -h		FLOAT	min p-value for exact test of HWE (two-tailed)
# -H	FLOAT	min p-value for exact test of excess heterozygotes
VARFILTER_OPTS="-v -a 0 -k $(( N_IND / 2 )) -u 2 -h 0 -H 1e-4"

# need a pileup file (not the best solution, could have used indexed bam or whatever)
samtools mpileup -DS $ALIGNS > $OUT/pileup

bcftools view $VARIANTS-raw.bcf | perl $NGSCLEAN/SNPcleaner.pl $VARFILTER_OPTS -A $ANC_SEQ -p $VARIANTS.filter_out1.vcf.bz2 -o $VARIANTS.TMP.vcf -G $OUT/pileup

# check depths for one chromosome only
cat $VARIANTS.TMP.vcf | grep chr1 | tr ";" "\t" | awk 'BEGIN{print "chr_pos\tdepth"} !/#/{sub("DP=","",$8); if(rand()<=0.05) print $1"_"$2"\t"$8}' > $OUT/chr1.sdepth

export R_LIBS_USER=~/brno2/sw/R-pkg211
Rscript --vanilla --slave $NGSCLEAN/get_depth_thres.R --in_file $OUT/chr1.sdepth --out_file $OUT/chr1.sdepth.fit.pdf --rnd_sample 20000 > $OUT/chr1.sdepth.depth_limits.tsv

time pbzip2 -p5 -dc /tmp/$ID.TMP.vcf.bz2 | tr ";" "\t" | awk 'BEGIN{print "chr_pos\tdepth"} !/#/{sub("DP=","",$8); if(rand()<=0.05) print $1"_"$2"\t"$8}' > $ID.sdepth
time Rscript --vanilla --slave $NGSCLEAN/get_depth_thresh.R --in_file $ID.sdepth --out_file $ID.sdepth.fit.pdf --rnd_sample 20000 > /tmp/$ID.depth_limits.tsv


#
# angsd
#
OUT=74-angsd
BAMS=$OUT/samples.txt
REF=51-liftover-all/lp2.fasta
SFSNAME=$OUT/folded

angsd -bam $BAMS -realSFS 1 -out $SFSNAME -anc $REF -GL 2 -fold 1 -nThreads $CPUS -minMapQ 1 -minQ 20
emOptim2 $SFSNAME.sfs $( cat $BAMS | wc -l ) -maxIter 100 -P $CPUS > $SFSNAME.sfs.em.ml

# per group
ls 63-smalt-all/ll*.bam > $OUT/ll.txt
ls 63-smalt-all/lm*.bam > $OUT/lm.txt
parallel -j 1 angsd -bam {} -realSFS 1 -out {.} -anc $REF -GL 2 -fold 1 -nThreads $CPUS -minMapQ 1 -minQ 20 ::: $OUT/l?.txt
parallel -j 1 'emOptim2 {.}.sfs $( cat {} | wc -l ) -maxIter 100 -P $CPUS > {.}.sfs.em.ml' ::: $OUT/l?.txt

# tajima's D
parallel -j 1 angsd -bam {} -GL 2 -out {.}-thetas -doThetas 1 -realSFS 1 -pest {.}.sfs.em.ml -fold 1 -anc $REF -nThreads $CPUS ::: $OUT/l?.txt

# major minor alele
parallel -j 1 angsd -bam {} -doMajorMinor 2 -doCounts 1 -out {.}-mm -nThreads $CPUS -minMapQ 1 -minQ 20 ::: $OUT/l?.txt

###
# use vcftools to calculate 'magic numbers'
# split species
$VCFLIB/vcfkeepsamples $OUT/variants-qual.vcf lu01 lu03 lu04 lu06 lu08 lu09 lu11 lu13 > $OUT/var-lm.vcf
$VCFLIB/vcfkeepsamples $OUT/variants-qual.vcf lu02 lu05 lu07 lu10 lu12 lu14 lu15 > $OUT/var-ll.vcf

vcftools --vcf $OUT/var-lm.vcf --out $OUT/var-lm --TsTv 1000 --window-pi 1000 --TajimaD 1000


# calculate N50 / N90 of the assembly
IN=22-newbler/454AllContigs.fna.filtered
QUANT=50
<$IN egrep -o 'length=[0-9]+' |
  sed 's/length=//' |
  mawk '{for(i=0;i<$0;i++) print $0;}' |
  sort -rn |
  perl -e '$d=.'$QUANT';@l=<>;print $l[int($d*$#l)]'


# fix terminal titles
# TODO: add MEM limit and session end time
echo -e "\[\033]0;$HOSTNAME: $PBS_NUM_PPN CPUS\007\]"

ANGSD=~/brno2/sw/angsd0.538
parallel "$ANGSD/misc/sfstools -sfsFile {} -nChr 15 -priorFile {}.sfs.ml -dumpBinary 1 > {}.sfs.norm" ::: $OUT/*.sfs

$ANGSD/misc/sfstools -sfsFile $OUT/ll.sfs -nChr 7 -priorFile $OUT/ll.sfs.em.ml -dumpBinary 1 > $OUT/ll.sfs.norm


###
#
# per site fst
tr " " "\n" > lm-samples <<< "lu01 lu03 lu04 lu06 lu08 lu09 lu11 lu13"
tr " " "\n" > ll-samples <<< "lu02 lu05 lu07 lu10 lu12 lu14 lu15"
vcftools --vcf 78-freebayes/variants-qual.vcf --out 78-freebayes/variants-qual --weir-fst-pop lm-samples --weir-fst-pop ll-samples

# use the mpileup variants, freebayes looks broken
IN=71-var-all/lp2-var-filtered.vcf.gz
vcftools --gzvcf $IN --out ${IN%%.*} --weir-fst-pop lm-samples --weir-fst-pop ll-samples

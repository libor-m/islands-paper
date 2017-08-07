library(tidyverse)
library(gtools)

# load data ####
load('data/daf_fst.RData')

#### QC of variant data ####

# does quality depend on read depth?
dvar %>%
  filter(qual < 999) %>%
  ggplot(aes(raw_depth, qual)) + geom_point(alpha=0.05) + scale_x_log10()

# qual ~ read depth plot
dvar %>%
  filter(qual < 999) %>%
  mutate(dpsum=ref_f + ref_r + alt_r + alt_f) %>%
  ggplot(aes(dpsum, qual)) + geom_point(alpha=0.1) + scale_x_log10()

# histogram of var qualities
dvar %>%
  filter(qual < 999) %>%
  ggplot(aes(qual)) + geom_histogram()

# 'gene' sizes in zebra finch
dvar %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(zf_len)) + geom_histogram()

# gene sizes per chromosome
dvar %>%
  group_by(chrom, contig_name) %>%
  summarize(zf_max=max(zf_max), zf_min=min(zf_min)) %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(chrom, zf_len)) + geom_boxplot() + coord_flip()

# max chrom sizes
dvar %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_max)) %>%
  ggplot(aes(chrom, chrom_len)) + geom_bar(stat="identity") + coord_flip()

# mvz annotations
dvar %>%
  ggplot(aes(mvz_seqpart)) +
  geom_histogram()

# filter out too wide contig targets
dvar %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(zf_len)) +
  geom_histogram()

# check whether 'chromosomes' are covered with variants
dfai <- read.delim("data-genome/lp2.fasta.fai",
                   header=F,
                   col.names=c("chrom", "len", "start", "x", "y"))

# get shared levels for final order
chrom_levs <- unique(c(levels(dvar$chrom), levels(dfai$chrom))) %>% mixedsort %>% rev
data.frame() %>%
  ggplot(aes(chrom)) +
  geom_bar(aes(y=len), stat="identity", data=dfai) +
  geom_point(aes(y=pos), colour="red", data=dvar) +
  coord_flip() +
  scale_x_discrete(limits=chrom_levs)
ggsave(file="results/var-coverage.png", width=200, height=290, units="mm")

#### QC of Fst values ####

# check check distribution of fst in badly mapped contigs
# this is not possible with current 'pre-filtered' approach
#daf %>%
#  ggplot(aes(fst, fill=zf_max - zf_min < 5e5)) +
#  geom_density(alpha=0.6, colour=NA)

daf %>%
  ggplot(aes(zf_max - zf_min < 5e5, fst)) +
  geom_boxplot() +
  coord_flip()
# maybe a bit lower fst values for the badly mapped exons

# check if the fst outside the mapped exons differ
daf %>%
  ggplot(daf, aes(fst, fill=is.na(zf_pos))) +
  geom_density(alpha=0.6, colour=NA)
# no apparent difference

# count not exactly mapped variants
daf %>%
  group_by(is.na(zf_pos)) %>%
  summarise(n())

# checking the zebra finch annotation, only 6 genes is longer than 500k (ensGenes)
daf %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len > 1e5) %>%
  .[,"contig_name"] %>%
  as.character %>%
  unique %>%
  length
# ->
# 1426 contigs are filtered out as 'too long' at 5e4
# 1005 at 1e5

# this filters out 16k additional variants
daf %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len > 5e4) %>%
  summarize(nvars=n())

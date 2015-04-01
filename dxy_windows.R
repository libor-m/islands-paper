# 
# look at dxy compared to Fst
#

setwd("c:/work/slavici-clanek/")

# granges must be loaded first beacuse it masks dplyr functions otherwise
# with some useless variants
library(GenomicRanges)

library(dplyr)
library(ggplot2)
library(gtools)

sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

# load the data and fix chromosome order
d <- read.delim("data/variant-table.tsv") %>% sortchrom
   
# load the dxy data
dxy <- read.delim('data/dxy.tsv', col.names=c("chrom", "pos", "ndiff", "ncomp", "dxy")) %>%
  filter(!is.na(ndiff))

# join on chrom, pos
da <- d %>% inner_join(dxy)

# pick 10 biggest chroms
bigchroms <- da %>% 
  filter(chrom != "chrUn") %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(pos)) %>%
  arrange(dplyr::desc(chrom_len)) %>%
  head(n=10) %>%
  .$chrom

# plot the variance ~ ncomparisons
dxy %>% ggplot(aes(ncomp, dxy)) + geom_point(alpha=0.1)
dxy %>% ggplot(aes(ncomp)) + geom_histogram()

# dot plot along chromosomes
# pick only reasonable values of ncomp based on previous plot
da %>% 
  filter(chrom %in% bigchroms, ncomp > 50) %>%
  ggplot(aes(pos, dxy)) + 
    geom_point(colour="#cccccc") + 
    facet_wrap(~chrom, ncol=1)

# use windows.R::find_variants to assign variants to sliding windows
daf_dxy <- da %>% filter(ncomp > 30)
ovr_dxy <- daf_dxy %>% find_variants

# calculate smoothed dxy
# improved since windows.R
smoothed_values <- function(hits, values) {
  hits %>%
    mutate(fst=values[subjectHits]) %>%
    group_by(chrom, zf_pos) %>%
    summarize(smooth=mean(fst, na.rm=T), nvars=n())
}

tdxy <- smoothed_values(ovr_dxy, daf_dxy$dxy) %>% dplyr::rename(dxy_smooth=smooth)
tdxy %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
    geom_line(colour="green") +
    facet_wrap(~chrom, ncol=1)
  
# get fst smoothers
dfst <- read.delim('data/lp2-var-filtered.weir.fst', col.names=c("chrom", "pos", "fst"))

daf_fst <- d %>% left_join(dfst) %>% filter(qual > 10)
ovr_fst <- daf_fst %>% find_variants
tfst <- smoothed_values(ovr_fst, daf_fst$fst) %>% dplyr::rename(fst_smooth=smooth)

# plot smooth fst and dxy
# looks i'm loosing the resolutino of fst by filtering the vars, go the other way around..
tdxy %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=dxy_smooth), colour="green") +
  geom_line(aes(y=fst_smooth), data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +
  facet_wrap(~chrom, ncol=1)

tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
  geom_line(colour="green") +
  facet_wrap(~chrom, ncol=1)

# weight dxy values by size of contig thei're in
# -> few big spikes probably due to missing data, let's check it
tdxy_scaled <- smoothed_values(ovr_dxy, daf_dxy$dxy / daf_dxy$contig_size) %>% dplyr::rename(dxy_smooth=smooth)
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=dxy_smooth), colour="green") +
  geom_line(aes(y=nvars), colour="#bbbbbb") +
  geom_line(aes(y=fst_smooth), data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +  
  facet_wrap(~chrom, ncol=1) +
  scale_y_log10() + 
  ggtitle("response of fst(blue) and dxy (green) to variant site density (gray)")

ggsave('results/metrics_var_density.pdf', width=20, height=16)

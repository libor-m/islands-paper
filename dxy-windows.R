# 
# look at dxy compared to Fst
#

setwd("c:/work/lab/slavici-clanek/")

# granges must be loaded first beacuse it masks dplyr functions otherwise
# with some useless variants
library(GenomicRanges)

library(dplyr)
library(ggplot2)
library(gtools)

sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

#### load data ####

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

#### sanity checks, overview ####

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

#### sliding window for Dxy ####

# use windows.R::find_variants to assign variants to sliding windows
daf_dxy <- da %>% filter(ncomp > 30)
ovr_dxy <- daf_dxy %>% find_variants

# evaluate the smoothed win
tdxy <- smoothed_values(ovr_dxy, daf_dxy$dxy) %>% rename(dxy_smooth=smooth)

# test plot
tdxy %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
    geom_line(colour="green") +
    facet_wrap(~chrom, ncol=1)
  
#### sliding window for Fst ####
dfst <- read.delim('data/lp2-var-filtered.weir.fst', col.names=c("chrom", "pos", "fst"))

daf_fst <- d %>% left_join(dfst) %>% filter(qual > 10)
ovr_fst <- daf_fst %>% find_variants
tfst <- smoothed_values(ovr_fst, daf_fst$fst) %>% rename(fst_smooth=smooth)

# plot smooth fst and dxy
# looks i'm loosing the resolution of fst by filtering the vars, go the other way around..
tdxy %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=dxy_smooth), colour="green") +
  geom_line(aes(y=fst_smooth), data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +
  facet_wrap(~chrom, ncol=1)

#### scale Dxy by variant density ####

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

tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
  geom_line(colour="green") +
  ylim(c(0, .0007)) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("per site dxy scaled by containing contig length")

ggsave('results/dxy_scaled.pdf', width=20, height=16)

#### rescale each variable to 0-1 scale, plot together (multitrack view) ####

# rescale numeric vector into (0, 1) interval
# clip everything outside the range 
rescale <- function(vec, lims=range(vec), clip=c(0, 1)) {
  # find the coeficients of transforming linear equation
  # that maps the lims range to (0, 1)
  slope <- (1 - 0) / (lims[2] - lims[1])
  intercept <- - slope * lims[1]
  
  xformed <- slope * vec + intercept
  
  # do the clipping
  xformed[xformed < 0] <- clip[1]
  xformed[xformed > 1] <- clip[2]
  
  xformed
}

# use brewer pallette for the overlapping lines
# to get same 'visual intensity'
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)), colour="Dxy")) +
  geom_line(aes(y=rescale(fst_smooth), colour="Fst"), data=tfst %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars), colour="nvars")) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1")

ggsave('results/metrics_var_density-scaled.pdf', width=20, height=16)


# check if the fst data is ok, find good ranges in histograms
# 
# fix the non-matching windows in the tboot and tfst .. ;(
# no time to calculate the bootstraps again
read.delim("data/tfst0.tsv") %>%
  select(-nvars) %>%
  inner_join(tfst) %>%
  rename(smooth=fst_smooth, boot=fst_boot) ->
  tfst_boot

read.delim("data/tdxy.tsv") %>%
  rename(smooth=dxy_smooth, boot=dxy_boot) ->
  tdxy_boot

# long formats for plotting
# (keep the wide ones for thresholding..)
library(tidyr)
tfst_boot %>%
  gather(type, fst, smooth, boot) ->
  tfst_boot_long

tdxy_boot %>%
  gather(type, dxy, smooth, boot) ->
  tdxy_boot_long

# histogram of fst and bootstraps
tfst_boot_long %>%
  ggplot(aes(fst, fill=type)) + 
  geom_histogram(alpha=0.6, position="identity")

# histogram of Dxy
tdxy_boot_long %>%
  ggplot(aes(dxy, fill=type)) + 
  geom_histogram(alpha=0.6, position="identity")

# one long table for all plotting - is this better?
t_long <- bind_rows(
  tfst_boot %>% gather(type, value, smooth, boot) %>% mutate(measure="Fst"),
  tdxy_boot %>% gather(type, value, smooth, boot) %>% mutate(measure="Dxy")
  )

# plot histograms for Fst and Dxy
t_long %>%
  ggplot(aes(value, fill=type)) +
  geom_histogram(position="identity", alpha=0.7) +
  facet_wrap(~measure)

# check nvars between measures
t_long %>%
  ggplot(aes(nvars, fill=measure)) +
  geom_histogram(position="identity", alpha=0.7)


tfst_boot_long %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, fst, colour=type)) + 
  geom_line() + 
  geom_point(y=0.2, colour="blue", shape=15, size=2, data=tfst_boot %>% filter(smooth > boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1)


# add bootstrap for fst
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)), colour="Dxy")) +
  geom_line(aes(y=rescale(smooth, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(fst_boot, lims=c(0, .25)), colour="Fst_boot"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars), colour="nvars")) +
  geom_point(y=1, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1")

ggsave('results/metrics_var_density-scaled.pdf', width=20, height=16)

# check if shifted tracks is more legible
# add means per value and chromosome
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)) + 1, colour="Dxy")) +
  geom_line(aes(y=rescale(dxy_boot, lims=c(0, .0007)) + 1, colour="Dxy_boot")) +
  geom_line(aes(y=rescale(smooth, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(fst_boot, lims=c(0, .25)), colour="Fst_boot"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars, lims=c(0, 700)) + 2, colour="nvars")) +
  geom_point(y=0, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .0007)) + 1, colour="Dxy"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(dxy_smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, 700)) + 2, colour="nvars"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=median(nvars))) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1") +
  theme(legend.title=element_blank())

ggsave('results/metrics_var_density-scaled-shift.pdf', width=20, height=16)


tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>% 
  group_by(chrom) %>% 
  summarize(cmean=mean(dxy_smooth) %>% rescale(c(0, .0007))) %>%
  ggplot(aes(chrom, cmean)) + geom_point() + geom_hline(aes(yintercept=0.3))

# TODO: center tracks around the mean?

# questions
# - how to deal with uneven sample coverage in contigs and across contigs
# - 


# calculate bootstrap for dxy
sample_dxy <- function(d)
  d %>% 
  select(chrom, zf_pos, dxy) %>%
  group_by(chrom) %>%
  mutate(dxy=sample(dxy)) %>%
  .[,"dxy"] %>%
  .[[1]]

reps <- 1:25000
lt <- sapply(reps, function(x) smoothed_values(ovr_dxy, sample_dxy(daf_dxy))$smooth)
save(lt, file='data/bootstraps-dxy.RData')
tdxy$dxy_boot <- apply(lt, 1, max)
tdxy_scaled$dxy_boot <- tdxy$dxy_boot
write.table(tdxy, "data/tdxy.tsv", sep="\t", quote=F, row.names=F)
rm(lt)
gc()


# and again, with dxy bootstrap
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)) + 1, colour="Dxy")) +
  geom_line(aes(y=rescale(smooth, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(fst_boot, lims=c(0, .25)), colour="Fst_boot"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars, lims=c(0, 700)) + 2, colour="nvars")) +
  geom_point(y=0, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .0007)) + 1, colour="Dxy"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(dxy_smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, 700)) + 2, colour="nvars"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=median(nvars))) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1") +
  theme(legend.title=element_blank())

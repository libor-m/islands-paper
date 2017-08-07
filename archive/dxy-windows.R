# 
# add Dxy to the measures
#

setwd("c:/work/lab/slavici-clanek/")

# granges must be loaded first beacuse it masks dplyr functions otherwise
# with some useless variants
library(GenomicRanges)

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)

# reverse version is good for sorting facets, because they're ordered bottom up usually
sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort))
sortchrom_r <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort %>% rev))

#### load data ####

# load the data and fix chromosome order
dvar <- read.delim("data/variant-table.tsv")
   
# load the dxy data
read.delim('data/dxy.tsv', 
           col.names=c("chrom", "pos", "ndiff", "ncomp", "dxy")) %>%
  filter(!is.na(ndiff)) ->
  ddxy

# join on chrom, pos 
# filter on some minimal variant quality
# and number of comparisons when calculating dxy
# (reflect number of sampled individuals for given site)
dvar %>%
  filter(zf_max - zf_min < 5e5) %>%
  filter(qual > 10) %>%
  inner_join(ddxy) %>%
  filter(ncomp > 30) ->
  daf_dxy

#### QC of Dxy data ####

# plot the variance ~ ncomparisons
ddxy %>% ggplot(aes(ncomp, dxy)) + geom_point(alpha=0.1)
ddxy %>% ggplot(aes(ncomp)) + geom_histogram()

# dot plot along chromosomes
# pick only reasonable values of ncomp based on previous plot
daf %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(pos, dxy)) + 
    geom_point(colour="#cccccc") + 
    facet_wrap(~chrom, ncol=1)

#### sliding window for Dxy ####

source('gr-bootstrap.R')

ovr_dxy <- daf_dxy %>% find_variants

# evaluate the smoothed win
tdxy <- smoothed_values(ovr_dxy, daf$dxy)

# test plot
tdxy %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, smooth)) + 
    geom_line(colour="green") +
    facet_wrap(~chrom, ncol=1)

# calculate the bootstrap
# calculate ~100 replicates and check the outcome
# before doing 25k
reps <- 1:100
lt <- sapply(reps, function(x) smoothed_values(ovr_dxy, rand_var(daf_dxy, "dxy"))$smooth)

smoothed_values(ovr_dxy, daf_dxy$dxy) %>% 
  mutate(boot = apply(lt, 1, max), measure="Dxy") ->
  tdxy_boot

tdxy_boot %>% 
  select(boot, smooth) %>% 
  gather(type, value) %>% 
  ggplot(aes(value, fill=type)) +
  geom_density(colour=NA, alpha=0.7)


# the full calculation
reps <- 1:25000
lt <- sapply(reps, function(x) smoothed_values(ovr_dxy, rand_var(daf_dxy, "dxy"))$smooth)
lt %>% save(file='data/bootstraps-dxy.RData')

# add bootstrap and measure tag to the real smoothed values
smoothed_values(ovr_dxy, daf_dxy$dxy) %>% 
  mutate(boot = apply(lt, 1, max), measure="Dxy") ->
  tdxy_boot

tdxy_boot %>% write.table("data/tdxy_boot.tsv", sep="\t", quote=F, row.names=F)

# clean up the huge table 
# not very useful in the end..
# https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14611
rm(lt)
gc()

#### sliding window for Fst ####

# load smooth and bootstrap produced in fst-windows.R
tfst_boot <- read.delim("data/tfst_boot.tsv")

# plot smooth fst and dxy
# looks i'm loosing the resolution of fst by filtering the vars, go the other way around..
tdxy %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, smooth)) + 
  geom_line(colour="green") +
  geom_line(data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +
  facet_wrap(~chrom, ncol=1)

#### scale Dxy by variant density ####

#FIXME: does the scaling make sense? shouldn't be bootstrap enough?

# weight dxy values by size of contig they're in
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

#### plot numbers together (multitrack view) ####

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

# long formats for plotting
# (keep the wide ones for thresholding..)
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

# version for tagged data sources
t_long <- bind_rows(
  tfst_boot %>% gather(type, value, smooth, boot),
  tdxy_boot %>% gather(type, value, smooth, boot)
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

# find quantiles
t_long %>% 
  group_by(paste(measure, type)) %>% 
  summarise(q01=quantile(value, 0.001), q99=quantile(value, 0.999))

# the multitrack plot with measures, bootstraps
# and 'islands'
# use 5 class diverging BrBG brewer palette
t_long %>%
  filter(chrom %in% bigchroms) %>%
  mutate(lineid=paste(measure, type),
         lineshift=measure %>% as.factor %>% as.numeric) %>%
  group_by(measure) %>%
  mutate(value=rescale(value, c(quantile(value, 0.001), quantile(value, 0.999)))) %>%
  ggplot(aes(zf_pos, value + lineshift, colour=lineid)) + 
  geom_line(aes(group=lineid)) + 
  geom_point(aes(colour="Fst smooth"), y=1,   shape=15, size=2, data=tfst_boot %>% filter(smooth > boot, chrom %in% bigchroms)) +
  geom_point(aes(colour="Dxy smooth"), y=0.8, shape=15, size=2, data=tdxy_boot %>% filter(smooth > boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  ylim(c(0.7, 3)) +
  scale_colour_manual(values=c(
    "Fst smooth"="#a6611a",
    "Fst"="#a6611a",
    "Fst boot"="#dfc27d",
    "Dxy smooth"="#018571",
    "Dxy"="#018571",
    "Dxy boot"="#80cdc1"
    ))

ggsave('results/multitrack.pdf', width=20, height=16)

#### spare parts (plotting) ####

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

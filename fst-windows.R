# 
# do some basic checks on the data and 
# calculate running window statistics for variants
#
setwd("c:/work/lab/slavici-clanek/")

library(dplyr)
library(ggplot2)
library(gtools)

#### load and merge data ####

# change chrom column in df to a factor with levels
# correctly sorted for plotting with ggplot
sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

dvar <- read.delim("data/variant-table.tsv")
dfst <- read.delim('data/lp2-var-filtered.weir.fst', col.names=c("chrom", "pos", "fst"))

# merge fst into variants information
# some rows in dfst are missing the value, filter those out
# filter out contigs badly mapped to zebra finch reference
# filter out low quality variants
dvar %>%
  inner_join(dfst) %>%
  filter(!is.na(fst)) %>%
  filter(zf_max - zf_min < 5e5) %>%
  filter(qual > 10) ->
  daf

# save filtered data
save(daf, file='data/daf.RData')
write.table(daf, file='data/daf.tsv', row.names=F, quote=F, sep="\t")

# load when resuming session
load('data/daf.RData')

#### QC of variant data ####

# start with some checks 
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

#### pick biggest chromosomes ####

# choose top n bigest chroms
dvar %>% 
  filter(chrom != "chrUn") %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_max)) %>%
  arrange(dplyr::desc(chrom_len)) %>%
  head(n=10) %>%
  .$chrom ->
  bigchroms

# alternatively, pick some..
bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")


#### load the bootstrap workers ####
source("gr-bootstrap.R")


#### timing checks ####
ptm <- proc.time()
ovr <- find_variants(daf)
proc.time() - ptm
# 11 seconds

ovr %>% View

ptm <- proc.time()
t <- smoothed_values(ovr, daf$fst)
proc.time() - ptm
# 0.5 seconds

ptm <- proc.time()
t <- smoothed_values(ovr, rand_var(daf, "fst"))
proc.time() - ptm
# <1 second
# that is 25k in few hours

# check if the result we got is reasonable
fst_plot <- function(x, chroms) 
  x %>% 
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, smooth)) +
  geom_line() +
  facet_wrap(~chrom, ncol=1)

t %>% fst_plot(bigchroms)

#### calculate the bootstrap ####
reps <- 1:25000
lt <- sapply(reps, function(x) smoothed_values(ovr, rand_fst(daf))$smooth)
lt %>% save(file='data/bootstraps-fst.RData')

# add bootstrap and measure tag to the real smoothed values
smoothed_values(ovr, daf$fst) %>% 
  mutate(boot = apply(lt, 1, max), measure="Fst") ->
  tfst_boot

tfst_boot %>% write.table("data/tfst.tsv", sep="\t", quote=F, row.names=F)

# clean up the huge table
rm(lt)
gc()

# a plot showing max bootstrapped value,
# the real smoothed value and the detected 'islands'
tm %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=fst_boot), colour="yellow") +
  geom_line(aes(y=fst_smooth), colour="blue") + 
  geom_point(aes(y=zf_pos), y=0.2, colour="blue", size=2, data=t %>% filter(fst_smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  ylim(c(-.1, .2)) +
  ggtitle("nightingale speciation islands as mapped to zebra finch chromosomes, 25k bootstrap")


# TODO - final filtering:
# - badly mapped contigs
# - low quality variants

# TODO - fst according to mvz seqpart

# - subset data - query only chromX data for chromX positions
# - parallelize from the outside

#### spare parts ####

# old slow and buggy approach
# do the means over 1 Mb window in zebra finch coordinates
# x is 'chrom', 'zf_pos'
# dd is the original dataset
window_average <- function(x, dd, window_half) {
  dd %>%
    filter(chrom == x[[1]]) %>%  # pick only variants on the same chromosome
    mutate(var_dist = abs(zf_pos - (x[[2]])) ) %>% # calculate distance to all variants from the current variant
    filter(var_dist < window_half) %>%
    summarize(mean_fst = mean(fst, na.rm=T)) %>% 
    .$mean_fst
}

# pick only zf positioned and real fst rows
daf <- da %>% filter(!is.na(zf_pos), !is.na(WEIR_AND_COCKERHAM_FST))
poslist <- split(daf[,c("chrom", "zf_pos")], seq_along(daf[,1]))

fst_smooth <- sapply(lpoints, window_average, daf, win_size/2)
daf$fst_smooth <-fst_smooth

# looks the calculation, apart from being hell slow,
# ended with many NAs in the results
chroms <- c("chr1", "chrZ")
daf %>%
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, WEIR_AND_COCKERHAM_FST)) + 
  geom_point(colour="#bbbbbb") + 
  geom_line(aes(y=fst_smooth, group=chrom)) +
  facet_wrap(~chrom, ncol=1)

daf %>%
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, fst_smooth)) + 
  geom_line() +
  facet_wrap(~chrom, ncol=1)

# move the data to python do do the window calculations
write.table(daf, file="data/vars-filtered-fst.tsv", sep="\t", row.names=F, quote=F)

daf <- read.delim('data/daf-fst.tsv') %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev)) 
ggplot(daf %>% filter(chrom == "chrZ"), aes(zf_pos, fst_1mb)) + geom_line()

#
# another failed approach
#
# a faster approach - an equal sampling along the chromosomes
# say 1k points along the longest chromosome
win_size <- 1e6
# required points
req_points <- 1e3
# ensure that the whole chrom is covered
stride <- as.integer(max(daf$zf_max) %>% { .  / max(. / win_size, req_points) })

# create equally spaced poins along all chromosomes, spaced by `stride`
dfpoints <- daf %>% 
  filter(chrom != "chrUn") %>%         # take only the known chromosomes
  mutate(chrom=as.character(chrom)) %>% # do not copy around the factor levels
  group_by(chrom) %>%                   # calculate chromosome sizes
  summarize(zf_size=max(zf_max)) %>%    # ..
  filter(zf_size > win_size) %>%        # pick chromosomes bigger than requested window
  split(seq_len(dim(.)[1])) %>%         # make data suitable for lapply
  lapply(function(x) 
    data.frame(chrom=x[[1]], zf_pos=seq(from=win_size/2, to=x[[2]], by=stride), stringsAsFactors=F)) %>%
  rbind_all                         # bind the dataframes together

# add the calculated values to the 'query' points
fst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, daf, win_size/2)
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth)
ggplot(dfst_smooth, aes(zf_pos, fst_smooth)) + geom_line() + facet_wrap(~chrom)

bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")
dfst_smooth %>% filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, fst_smooth)) + geom_line() + facet_wrap(~chrom, ncol=1)

# qual - fst plot
daf %>% filter(chrom %in% bigchroms, qual < 900) %>%
  ggplot(aes(qual, fst)) + geom_point(alpha=0.1) + facet_wrap(~chrom)

# point + line plot
# does not make sens in ggplot, because there is no dual axes in ggplot, and the
# scales are different

# fst smooth on quality filtered data
dafq <- daf %>% filter(qual > 10)
fst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, dafq, win_size/2)
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth)

# randomize the fst values in each chromosome
dfr <- dafq %>% 
  select(chrom, zf_pos, fst) %>%
  group_by(chrom) %>%
  mutate(fst=sample(fst))

# 
ptm <- proc.time()
rfst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, dfr, win_size/2)
proc.time() - ptm
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth, fst_random = as.numeric(rfst_smooth))

dfst_smooth %>% filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=fst_smooth), colour="blue") + 
  geom_line(aes(y=fst_random), colour="yellow") +
  facet_wrap(~chrom, ncol=1)

# 
# do some basic checks on the data and 
# calculate running window statistics for variants
#
setwd("c:/work/slavici-clanek/")

# load the data and fix chromosome order
library(gtools)
d <- read.delim("data/variant-table.tsv") %>% 
  mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev)) 

library(dplyr)
library(ggplot2)

# start with some checks 
# does quality depend on read depth?
d %>% 
  filter(qual < 999) %>%
  ggplot(aes(raw_depth, qual)) + geom_point(alpha=0.05) + scale_x_log10()

# qual ~ read depth plot
d %>% 
  filter(qual < 999) %>%
  mutate(dpsum=ref_f + ref_r + alt_r + alt_f) %>%
  ggplot(aes(dpsum, qual)) + geom_point(alpha=0.1) + scale_x_log10()

# histogram of var qualities
d %>% 
  filter(qual < 999) %>%
  ggplot(aes(qual)) + geom_histogram()

# 'gene' sizes in zebra finch
d %>% 
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(zf_len)) + geom_histogram()

# gene sizes per chromosome
d %>%
  group_by(chrom, contig_name) %>%
  summarize(zf_max=max(zf_max), zf_min=min(zf_min)) %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%  
  ggplot(aes(chrom, zf_len)) + geom_boxplot() + coord_flip()

# max chrom sizes
d %>% 
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_max)) %>%
  ggplot(aes(chrom, chrom_len)) + geom_bar(stat="identity") + coord_flip()

# mvz annotations
d %>%
  ggplot(aes(mvz_seqpart)) + geom_histogram()


# merge in fst data
dfst <- read.delim('data/lp2-var-filtered.weir.fst')
da <- d %>% left_join(dfst, c("chrom" = "CHROM", "pos" = "POS"))

# save for further reuse
colnames(da)[19] <- "fst"
# get rid of positions without fst
daf <- da %>% filter(!is.na(fst))

save(daf, file='data/daf.RData')
write.table(daf, file='data/daf.tsv', row.names=F, quote=F, sep="\t")

# check if the fst outside the mapped exons differ
ggplot(daf, aes(fst, fill=is.na(zf_pos))) + geom_density(alpha=0.6, colour=NA)

# plot few chromosomes
chroms <- c("chr1", "chrZ")
da %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST), chrom %in% chroms) %>%
  ggplot(aes(pos, WEIR_AND_COCKERHAM_FST)) + geom_point(colour="#bbbbbb") + facet_wrap(~chrom, ncol=1)

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


# attempts to speed things up (naive approach takes 180 s on single core)
# - range query
library(IRanges)
ir <- IRanges(start=ifelse(is.na(daf$zf_pos), daf$zf_min, daf$zf_pos), 
              end=ifelse(is.na(daf$zf_pos), daf$zf_max, daf$zf_pos + 1))
rd <- RangedData(ir, space=daf$chrom, daf$fst )
library(GenomicRanges)

gr <- GRanges(seqnames=daf$chrom, ranges=ir, mcols=daf[,"fst"])
gf <- GIntervalTree(gr)

irq <- IRanges(start=dfpoints$zf_pos - win_size / 2, end=dfpoints$zf_pos + win_size / 2)
gq <- GRanges(seqnames=dfpoints$chrom, ranges=irq)

ovr <- findOverlaps(gq, gf)
ovr %>% as.data.frame %>% arrange(queryHits) %>% View

# 
t <- ovr %>% as.data.frame %>% 
  mutate(chrom=as.factor(seqnames(gf)[subjectHits]), 
         fst=mcols(gf)[subjectHits,], 
         zf_pos=gq %>% ranges %>% start %>% {.[queryHits] + win_size / 2}  ) %>%
  group_by(chrom, zf_pos) %>%
  summarize(fst_smooth=mean(fst))

# this works, and is faster..
# wrap to functions

# - subset data - query only chromX data for chromX positions
# - parallelize from the outside


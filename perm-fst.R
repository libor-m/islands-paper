#
# detect islands with a permutation based model
#

library(tidyverse)
library(gtools)

sortchrom <- function(df) df %>% ungroup %>%mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

#### load and merge data ####

# load
load('data/daf_fst.RData')

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
source("lib/bootstrap.R")

# a few checks

system.time(ovr <- find_variants(daf_fst))
# 11 seconds
# .37 seconds with data.table

ovr %>% View

system.time(
  t <- smoothed_values(ovr, daf$fst))
# 0.5 seconds

system.time(
  t <- smoothed_values(ovr, rand_var(daf, "fst")))
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
# group the variants
daf %>%
  mutate(AZ = grepl("^chrZ", chrom)) %>%
  group_by(AZ) ->
  daf_az

# reps <- 1:250
reps <- 1:25000
system.time(
  lt <- sapply(reps, function(x) smoothed_values(ovr, rand_var(daf_az, "fst"))$smooth))
#    user  system elapsed
# 4951.27  253.64 5247.94

# save the precious result immediately ;)
lt %>% save(file = 'data/bootstraps-fst.RData')

# add bootstrap and measure tag to the real smoothed values
smoothed_values(ovr, daf$fst) %>%
  mutate(boot = apply(lt, 1, max),
         boot_q99 = apply(lt, 1, quantile, .99, type = 8),
         boot_q95 = apply(lt, 1, quantile, .95, type = 8),
         measure = "Fst") ->
  tfst_boot

tfst_boot %>% write.table("data/tfst_boot.tsv", sep="\t", quote=F, row.names=F)

# clean up the huge table
rm(lt)
gc()

# a plot showing max bootstrapped value,
# the real smoothed value and the detected 'islands'
tfst_boot %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) +
  geom_line(aes(y=boot), colour="yellow") +
  geom_line(aes(y=smooth), colour="blue") +
  geom_point(aes(y=zf_pos), y=0.2, colour="blue", size=2, data=tfst_boot %>% filter(smooth > boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1, ) +
  ylim(c(-.1, .2)) +
  ggtitle("nightingale speciation islands as mapped to zebra finch chromosomes, 25k bootstrap")


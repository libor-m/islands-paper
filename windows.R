# calculate running window statistics for variants
d <- read.delim("data/variant-table.tsv")

library(dplyr)
library(ggplot2)

# start with some checks 
# does quality depend on read depth?
d %>% 
  filter(qual < 999) %>%
  ggplot(aes(raw_depth, qual)) + geom_point()

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

# mvz annotations
d %>%
  ggplot(aes(mvz_seqpart)) + geom_histogram()



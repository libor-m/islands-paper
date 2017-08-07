library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

#### load the data ####

read_delim('data/d_wins.tsv', delim=" ") ->
  dwins

# dxy in lp coordinates
read_tsv('data/dxy.tsv',
         col_names=c("chrom", "pos", "ndiff", "ncomp", "dxy")) %>%
  filter(!is.na(ndiff)) ->
  ddxy

# read all other info on the variants
read_tsv("data/variant-table.tsv") -> dvar

# join on chrom, pos
# filter on some minimal variant quality
# and number of comparisons when calculating dxy
# (reflect number of sampled individuals for given site)
dvar %>%
  filter(zf_max - zf_min < 5e5) %>%
  filter(qual > 10) %>%
  inner_join(ddxy, by=c("chrom", "pos")) %>%
  filter(ncomp > 30) %>%
  mutate(zf_pos = ifelse(is.na(zf_pos),
                         (zf_max - zf_min) / 2,
                         zf_pos) %>% as.integer) ->
  daf_dxy

#### separate the variants according to in/out of Fst island ####

source('lib/interval-tools.R')

# convert `zf_pos` - position in zebra finch reference to intervals
daf_dxy %>%
  mutate(start = zf_pos, end = zf_pos) ->
  ddxy_int

# dt.intersect(ddxy_int, dwins %>% filter(measure == "Fst")) ->  dint
# there is some problem with setkey() in the intersection
# let's do it here
dtxy <- as.data.table(ddxy_int)
setkey(dtxy, chrom, start, end)

dtwin <- as.data.table(dwins %>% filter(measure == "Fst"))
setkey(dtwin, chrom, start, end)

# pick inside variants
dt.intersect(dtxy, dtwin) ->  d_in_islands

# create a fake genome
# with chromosomes ending after the last variant
dvar %>%
  group_by(chrom) %>%
  summarise(len = max(zf_max) + 1) ->
  genome

dwins %>%
  filter(measure == "Fst") %>%
  interval.complement(genome) ->
  fst_non_windows

# check if the complement is working
bind_rows(fst_non_windows %>% mutate(measure = "non-window"),
          dtwin %>% mutate(measure = "Fst window")) %>%
  mutate(id = row_number()) %>%
  gather(type, pos, start, end) %>%
  ggplot(aes(pos, colour=measure)) +
  geom_line(aes(y=measure, group=id), size=3) +
  facet_wrap(~chrom, ncol = 2, switch="y")
ggsave('results/complement-check.pdf', width=9, height=12)

dtnonwin <- as.data.table(fst_non_windows)
setkey(dtnonwin, chrom, start, end)

dt.intersect(dtxy, dtnonwin) -> d_out_islands

#### now put boths sets together and compare ####
bind_rows(d_in_islands %>% select(chrom, pos = start, dxy) %>% mutate(group="IN"),
          d_out_islands %>% select(chrom, pos = start, dxy) %>% mutate(group="OUT")) ->
  dxy_islands

dxy_islands %>%
  ggplot(aes(group, dxy, fill=group)) +
  geom_boxplot()
ggsave('results/dxy-inout-boxplot.pdf', width=9, height=12)
ggsave('results/dxy-inout-boxplot.png', width=9, height=12, dpi=72)

dxy_islands %>%
  ggplot(aes(dxy, fill=group)) +
  geom_density(colour = NA, alpha = 0.7)
ggsave('results/dxy-inout-density.pdf', width=12, height=9)
ggsave('results/dxy-inout-density.png', width=12, height=9, dpi=72)

# significance testing
t.test(d_in_islands$dxy, d_out_islands$dxy)
# t = 4.6712, df = 3314.368, p-value = 3.113e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01284956 0.03143910
# sample estimates:
#   mean of x mean of y
# 0.2995938 0.2774495


# compare the overlaps of the Fst and Dxy windows
# simply by taking intersection length / union length
bind_rows(dwins %>%
            filter(measure != "both") %>%
            reduce_intervals %>%
            mutate(measure = "union"),
          dwins %>%
            filter(measure == "both") %>%
            mutate(measure = "intersection")) %>%
  group_by(measure) %>%
  mutate(len=end - start) %>%
  summarise(len = sum(len))

# calculate the ratio (Jaccard index)
84456478 / 303645680
# 0.2781415

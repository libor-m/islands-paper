#
# calculate some basic stats on the data
# to report in the paper
#
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

read_delim("data-genome/vars-melted.tsv",
           delim = " ",
           col_names = c("chrom", "pos", "individual", "allele", "seq")) ->
  d

read_delim("data/populations",
           delim = " ",
           col_names = c("individual", "species")) ->
  ds

d %>%
  left_join(ds, by="individual") ->
  dv

dv %>% write.table("data-genome/vars-melted-species.tsv",
                   sep = "t",
                   quote = F,
                   row.names = F,
                   col.names = F)


# histogram of coverage calculated with sort | uniq
read_tsv("data/alldup.depth.hist", col_names = c("count", "coverage")) ->
  dh

dh %>% 
  ggplot(aes(coverage, count)) + 
  geom_point() + 
  scale_x_log10()
ggsave('results/lp2-coverage.png', width=10, height=10, dpi=72)

read_tsv("data-genome/lp2-var-filtered.depths", col_names = c("depth")) ->
  dvh

dvh$depth %>% mean

dvh %>%
  ggplot(aes(depth)) +
  geom_histogram(bins=60) + 
  scale_x_log10()
ggsave('results/lp2-var-filtered-DP.png', width=6, height=6, dpi=72)

# how many islands
read_tsv("data/d_wins.tsv") ->
  dislands

dislands %>%
  group_by(measure) %>%
  summarise(n = n(), size = sum(end - start))

# the window size is a bit dominated by the 1M window around each 
# positive value, but for Fst it looks 
dislands %>%
  filter(measure != "both") %>%
  mutate(length=end - start) %>%
  ggplot(aes(length)) + 
  geom_histogram() +
  facet_wrap(~measure)

# stats on the islands
dislands %>%
  # filter(measure != "both") %>%
  mutate(length=end - start) %>%
  group_by(measure) %>%
  summarise(count = n(), 
            covered = sum(length),
            mean = mean(length),
            min=min(length),
            max=max(length)) %>%
  write.table("clipboard", sep="\t", row.names=F)


# zf chromosome length ----

daf_fst %>%
  group_by(autosome = !(chrom %in% c("chrZ", "chrZ_random"))) %>%
  summarize(len = max(zf_max))

#  autosome       len
# 1    FALSE  72831270
# 2     TRUE 175120622
#

# check if there is more bases (fractional) covered in autosomes or chrZ
dislands %>%
  group_by(measure, autosome = !(chrom %in% c("chrZ", "chrZ_random"))) %>%
  summarize(len = sum(end - start)) %>%
  mutate(autosome = ifelse(autosome, "autosome", "chrZ")) %>%
  spread(autosome, len) %>%
  mutate(autosome = autosome / 175120622,
         chrZ = chrZ / 72831270)

# chrom size versus bases covered in islands

daf_fst %>%
  group_by(chrom) %>%
  summarize(chrom_size = max(zf_max)) %>%
  left_join(dislands %>% 
              filter(measure == "Fst") %>%
              group_by(chrom) %>%
              summarise(in_fst_wins = sum(end - start)),
            by = "chrom") %>%
  filter(!is.na(in_fst_wins)) %>%
  ggplot(aes(chrom_size, in_fst_wins)) +
  geom_text(aes(label = chrom)) +
  geom_smooth(method = "lm")

# count fst windows
daf_fst %>%
  group_by(chrom) %>%
  summarize(chrom_size = max(zf_max)) %>%
  left_join(dislands %>% 
              filter(measure == "Fst") %>%
              group_by(chrom) %>%
              summarise(fst_wins = n()),
            by = "chrom") %>%
  filter(!is.na(fst_wins)) %>%
  ggplot(aes(chrom_size, fst_wins)) +
  geom_text(aes(label = chrom)) +
  geom_smooth(method = "lm")

# TODO: chisq

# average Fst
daf_fst %>%
  group_by(autosome = !(chrom %in% c("chrZ", "chrZ_random"))) %>%
  summarise(meanFst = mean(fst), medFst = median(fst))

daf_fst %>%
  mutate(autosome = !(chrom %in% c("chrZ", "chrZ_random"))) %>%
  ggplot(aes(fst, fill=autosome)) + 
  geom_density(colour=NA, alpha=0.7, adjust=3)

  

#
# auxiliary plots for publication/presentation
#

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)

#### load data ####
tfst_boot <- read.delim("data/tfst_boot.tsv")
tdxy_boot <- read.delim("data/tdxy_boot.tsv")

t_boot <- bind_rows(tfst_boot, tdxy_boot)

#### helpers ####

sortchrom <- function(df) df %>% ungroup %>%mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

t_boot %>% 
  filter(chrom != "chrUn") %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_pos)) %>%
  arrange(dplyr::desc(chrom_len)) %>%
  head(n=10) %>%
  .$chrom ->
  bigchroms

#### dxy and fst for bigger chromosomes ####

t_boot %>%
  filter(chrom %in% bigchroms) %>%
  mutate(chrom_type = ifelse(chrom == "chrZ", "sex chromosome", "autosome")) %>%
  ggplot(aes(chrom, smooth)) +
  geom_boxplot(aes(fill=chrom_type)) +
  facet_wrap(~ measure, ncol = 1, scales="free_y") +
  ylab("")
ggsave("results/fst.dxy-boxplot-top10.pdf", width=200, height=150, units="mm")

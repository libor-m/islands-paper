#### figure 1 for the paper ####
# show ratio of real fst vs bootstrap, so the peaks are
# only in the areas of high differentiation

library(tidyverse)
library(gtools)

# reverse version is good for sorting facets, because they're ordered bottom up usually
sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort))
sortchrom_r <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort %>% rev))

bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")

read_tsv("data/tfst_boot.tsv") %>%
  select(chrom, zf_pos, Fst=smooth, Bootstrap=boot_q99) ->
  tfst_boot

# get ddxy and dovr from dxy-in-fst-islands.R

ddxy$dxy_rel %>% mean -> dxy_genomic_mean

dovr %>%
  group_by(chrom, start, end, id) %>%
  summarise(dxy_rel = mean(dxy_rel)) %>%
  gather(type, zf_pos, start, end) %>%
  ungroup ->
  dovr_long

tfst_boot %>%
  filter(chrom %in% bigchroms) %>%
  sortchrom_r %>%
  ggplot(aes(zf_pos)) +
  geom_line(aes(y = Fst / Bootstrap), colour = "#999999") +
  geom_hline(yintercept = 1, colour="firebrick", alpha=0.6) +
  geom_line(aes(group=id, colour=dxy_rel > dxy_genomic_mean),
            y = -.1,
            size = 3,
            data = dovr_long %>% filter(chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol = 1, strip.position = "left")
  # facet_grid(cgroup ~ chrom, space = "free_x", scales = "free_x") +
  # guides(colour = FALSE)

# do all chromosomes and rearrange in illustrator
tfst_boot %>%
  filter(!grepl("random", chrom)) %>%
  sortchrom %>%
  ggplot(aes(zf_pos)) +
  geom_line(aes(y = Fst / Bootstrap),
            colour = "#999999") +
  geom_hline(yintercept = 1,
             colour = "firebrick",
             alpha = 0.6) +
  geom_line(aes(group = id,
                colour = dxy_rel > dxy_genomic_mean),
            y = -.1,
            size = 3,
            data = dovr_long %>%
                     filter(!grepl("random", chrom)) %>%
                     sortchrom) +
  facet_wrap(~chrom,
             ncol = 1,
             strip.position = "left") +
  xlim(0, 3e8) +
  scale_colour_manual(values = c("FALSE" = "#aaaaaa",
                                 "TRUE" = "firebrick"),
                      guide = "none")

# 31 panels, we want ~3 cm per panel, that is 31 * 30
ggsave('results/figure1-src-rev1.pdf', width = 290, height = 900, unit = "mm")

tfst_boot %>%
  filter(!grepl("random|LG", chrom)) %>%
  sortchrom %>%
  ggplot(aes(zf_pos)) +
  geom_line(aes(y = Fst),
            colour = "firebrick") +
  geom_line(aes(y = Bootstrap),
            colour = "#999999") +
  geom_line(aes(group = id,
                colour = dxy_rel > dxy_genomic_mean),
            y = -.1,
            size = 3,
            data = dovr_long %>%
                     filter(!grepl("random", chrom)) %>%
                     sortchrom) +
  facet_wrap(~chrom,
             ncol = 1,
             strip.position = "left") +
  xlim(0, 3e8) +
  ylim(-0.1, 0.3) +
  scale_colour_manual(values = c("FALSE" = "#aaaaaa",
                                 "TRUE" = "firebrick"),
                      guide = "none")

ggsave('results/figure1-src-rev1-both.pdf', width = 290, height = 900, unit = "mm")

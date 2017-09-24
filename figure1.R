#### figure 1 for the paper ####
# show ratio of real fst vs bootstrap, so the peaks are
# only in the areas of high differentiation

library(tidyverse)
library(gtools)

# local library ----
# reverse version is good for sorting facets, because they're ordered bottom up usually
sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort))
sortchrom_r <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% unique %>% mixedsort %>% rev))

# setup data ----
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

# plot Fst in single line ----

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

# plot Fst ans Bootstrap separately ----

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

# wrap chromosomes programatically ----

# the outcome is a fixed width plot, where chromosomes 
# are ordered by mixed sort and start at round positions, 
# namely 300 MB wide, chromosomes starting each 50 MB
tfst_boot %>%
  count(chrom) %>%
  filter(n >= 30,
         !grepl("random|LG", chrom)) %>%
  .$chrom ->
  display_chroms

tfst_boot %>%
  filter(chrom %in% display_chroms) %>%
  sortchrom %>%
  arrange(chrom, zf_pos) ->
  tfst_plotdata

block_size <- 5e7
panel_size <- 3.5e8

tfst_plotdata %>%
  # each chrom occupies space sized in multiples of block_size
  group_by(chrom) %>%
  summarise(size = zf_pos %>% max,
            round_size = size %>% plyr::round_any(block_size, ceiling)) %>%
  # each chromosome has to fit on single row, move it to next panel otherwise
  {reduce2(.$chrom, .$round_size, function(acc, chrom, size) {
    next_break <- plyr::round_any(acc$tail, panel_size, ceiling)
    start <- if_else(acc$tail + size > next_break, next_break, acc$tail)
    list(tail = start + size,
         starts = bind_rows(
           acc$starts,
           list(chrom = chrom, start = start, round_size = size)
         )
    )
  }, .init = list(tail = 0, starts = list()))} %>%
  .$starts %>%
  mutate(
    offset_total = start,
    end_total = offset_total + round_size,
    panel = (end_total - 1) %/% panel_size,
    offset_panel = pmax(0, offset_total - (panel * panel_size))) %>% 
  select(chrom, panel, offset_panel) ->
  chrom_panels

# Fst, measure and q99 line ----
tfst_plotdata %>%
  left_join(chrom_panels, by = "chrom") %>%
  mutate(pos = zf_pos + offset_panel) %>%
  ggplot(aes(pos, group = chrom)) +
  geom_line(aes(y = Fst),
            colour = "firebrick") +
  geom_line(aes(y = Bootstrap),
            colour = "#999999") +
  
  # islands marks below
  geom_line(aes(group = id,
                colour = dxy_rel > dxy_genomic_mean),
            y = -.1,
            size = 3,
            data = dovr_long %>%
                     filter(chrom %in% display_chroms) %>%
                     left_join(chrom_panels, by = "chrom") %>%
                     mutate(pos = zf_pos + offset_panel)) +
  
  # chromosome labels with hacky width of the rect
  geom_rect(aes(x=offset_panel, xmin = offset_panel - 5e6, xmax = offset_panel),
            ymin = -Inf,
            ymax = Inf,
            fill = "#aaaaaa",
            data = chrom_panels) +
  geom_text(aes(offset_panel, label = chrom),
            angle = 90,
            y = 0,
            hjust = "left",
            vjust = -0.3,
            data = chrom_panels) +
  
  facet_wrap(~panel,
             ncol = 1,
             strip.position = "right") +
  xlim(-5e6, panel_size) +
  ylim(-0.1, 0.3) +
  scale_colour_manual(values = c("FALSE" = "#aaaaaa",
                                 "TRUE" = "firebrick"),
                      guide = "none")
  
ggsave('results/figure1-src-rev1-both.pdf', width = 290, height = 210, unit = "mm")

# Fst ratio ----
tfst_plotdata %>%
  left_join(chrom_panels, by = "chrom") %>%
  mutate(pos = zf_pos + offset_panel) %>%
  ggplot(aes(pos, group = chrom)) +
  
  # the ratio line
  # coloring by fst does not look nice (would need additional 
  # points at y = 1 to end the coloring..)
  geom_line(aes(y = Fst / Bootstrap),
            colour = "#999999") +
  
  geom_hline(yintercept = 1,
             colour = "firebrick",
             alpha = 0.6) +
  
  # islands marks below
  geom_line(aes(group = id,
                colour = dxy_rel > dxy_genomic_mean),
            y = -.1,
            size = 3,
            data = dovr_long %>%
                     filter(chrom %in% display_chroms) %>%
                     left_join(chrom_panels, by = "chrom") %>%
                     mutate(pos = zf_pos + offset_panel)) +
  
  # chromosome labels with hacky width of the rect
  geom_rect(aes(x=offset_panel, xmin = offset_panel - 5e6, xmax = offset_panel),
            ymin = -Inf,
            ymax = Inf,
            fill = "#aaaaaa",
            data = chrom_panels) +
  geom_text(aes(offset_panel, label = chrom),
            angle = 90,
            y = 0,
            hjust = "left",
            vjust = -0.3,
            data = chrom_panels) +
  
  facet_wrap(~panel,
             ncol = 1,
             strip.position = "right") +
  xlim(-5e6, panel_size) +
  scale_colour_manual(values = c("FALSE" = "#aaaaaa",
                                 "TRUE" = "firebrick"),
                      guide = "none")
  
ggsave('results/figure1-src-rev1-ratio.pdf', width = 290, height = 210, unit = "mm")

# dxy both values ----
read_tsv("data/tdxy_boot.tsv") %>%
  filter(measure == "Dxy") %>%
  select(chrom, zf_pos, Dxy = smooth, Bootstrap = boot_q99) ->
  tdxy_boot

tdxy_boot %>%
  filter(chrom %in% display_chroms) %>%
  left_join(chrom_panels, by = "chrom") %>%
  mutate(pos = zf_pos + offset_panel) %>%
  ggplot(aes(pos, group = chrom)) +
  geom_line(aes(y = Dxy),
            colour = "seagreen") +
  geom_line(aes(y = Bootstrap),
            colour = "#999999") +
  
  # islands marks below
  geom_line(aes(group = id,
                colour = dxy_rel > dxy_genomic_mean),
            y = 0,
            size = 3,
            data = dovr_long %>%
                     filter(chrom %in% display_chroms) %>%
                     left_join(chrom_panels, by = "chrom") %>%
                     mutate(pos = zf_pos + offset_panel)) +
  
  # chromosome labels with hacky width of the rect
  geom_rect(aes(x=offset_panel, xmin = offset_panel - 5e6, xmax = offset_panel),
            ymin = -Inf,
            ymax = Inf,
            fill = "#aaaaaa",
            data = chrom_panels) +
  geom_text(aes(offset_panel, label = chrom),
            angle = 90,
            y = 0.003,
            hjust = "left",
            vjust = -0.3,
            data = chrom_panels) +
  
  facet_wrap(~panel,
             ncol = 1,
             strip.position = "right") +
  xlim(-5e6, panel_size) +
  ylim(0, 0.01) +
  scale_colour_manual(values = c("FALSE" = "#aaaaaa",
                                 "TRUE" = "firebrick"),
                      guide = "none")
  
ggsave('results/figure1-src-rev1-dxy.pdf', width = 290, height = 210, unit = "mm")

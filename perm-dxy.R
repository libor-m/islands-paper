#
# detect islands with a permutation based model
# - it's possible to use the per gene dxy, because 
#   contigs are smaller than the windows, thus can be aggregated
#   in the same way as single point mutations
#
library(tidyverse)

#### load data ####
# setwd("c:/work/lab/slavici-clanek/")

read_tsv('data-gene-dxy/contigs_dxy_avg_to_fst_isl.txt') -> ddxy

# make the new data by vasek as similar as possible to the previous ones
read_tsv('data-gene-dxy/dxy_abs_by_contig.tab') %>%
  rename(
    zf_contig_chrom = ZF_CONTIG_CHROM,
    zf_contig_start = ZF_CONTIG_START,
    zf_contig_end = ZF_CONTIG_END,
    contig_id = CONTIG,
    contig_length = SIZE,
    dxy_abs = DXY_ABS,
    no_overlapping_islands = FST_ISLAND) ->
  ddxy

#### QC and visual analytics ####

# gene length seems to be bounded to 50k
ddxy %>%
  ggplot(aes(zf_contig_chrom, zf_contig_end - zf_contig_start)) +
  geom_boxplot()

# contig length vs gene length
ddxy %>%
  ggplot(aes(contig_length, zf_contig_end - zf_contig_start)) +
  geom_point()

# dxy vs contig length
# looks like there is a non-linear dependency, dxy gets
# saturated at contig len ~1k
ddxy %>%
  ggplot(aes(contig_length, dxy_abs)) +
  geom_point(alpha = 0.1) +
  scale_y_log10() +
  xlim(0, 2500)
ggsave(
  'results/dxy-contig-length.png',
  width = 8,
  height = 10,
  dpi = 72
)

ddxy %>%
  ggplot(aes(contig_length)) +
  geom_histogram()

#### the bootstrap ####

source('lib/bootstrap.R')

# rename columns to work with find_variants
# and run find_variants
# `dxy_rel` should be calculated by dividing with 
#   the amount of bases with the same filtering criteria as the SNPs,
#   but the ratio of filtered bases to all bases will be ~constant, 
#   that is `dxy_rel` will differ only by a constant
#
ddxy %>%
  dplyr::rename(
    chrom = zf_contig_chrom) %>%
  mutate(
    zf_min = zf_contig_start %>% as.numeric,
    zf_max = zf_contig_end %>% as.numeric,
    zf_pos = zf_min + ((zf_max - zf_min) / 2),
    dxy_rel = dxy_abs / contig_length ) %>%
  select(
    -zf_contig_start,
    -zf_contig_end) ->
  ddxy_fix

ddxy_fix %>%
  find_variants ->
  ovr_dxy

# evaluate the smoothed win
tdxy <- smoothed_values(ovr_dxy, ddxy$dxy_abs)

# test plot
bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")

tdxy %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, smooth)) +
  geom_line(data = tdxy_orig %>% filter(chrom %in% bigchroms), colour = "gray70") +
  geom_line(colour = "green") +
  facet_wrap(~chrom, ncol = 1, strip.position = "left")
ggsave("results/dxy-gene-old_new.pdf", width = 9, height = 7)

ddxy_fix %>%
  mutate(AZ = grepl("^chrZ", chrom)) %>%
  group_by(AZ) ->
  ddxy_az

# calculate the bootstrap
# calculate ~100 replicates and check the outcome
# before doing 25k
reps <- 1:100
reps <- 1:25000

system.time(
  lt <- sapply(reps, function(x) smoothed_values(ovr_dxy, rand_var(ddxy_az, "dxy_rel"))$smooth))
#     user   system  elapsed
# 12792.97    50.87 12848.47

save(lt, file = 'data-gene-dxy/bootstraps-dxy-gene-rel.RData')
# load bootstrap data if ever needed again;)
load('data-gene-dxy/bootstraps-dxy-gene-rel.RData')

# relative, quantiles
smoothed_values(ovr_dxy, ddxy_fix$dxy_rel) %>%
  mutate(boot_max = apply(lt, 1, max),
         boot_q75 = apply(lt, 1, quantile, probs = 0.75),
         boot_q95 = apply(lt, 1, quantile, probs = 0.95),
         boot_q99 = apply(lt, 1, quantile, probs = 0.99),
         measure = "Dxy") ->
  tdxy_boot

tdxy_boot %>%
  write.table('data-gene-dxy/dxy-gene-bootstrap.tsv',
              sep = "\t",
              quote = F,
              row.names = F)
tdxy_boot %>%
  write.table('data/tdxy_boot.tsv',
              sep = "\t",
              quote = F,
              row.names = F)

read_tsv('data-gene-dxy/dxy-gene-bootstrap.tsv') -> tdxy_boot

tdxy_boot %>%
  gather(type, value, smooth, boot_max:boot_q99) %>%
  ggplot(aes(value, fill=type)) +
  geom_density(colour=NA, alpha=0.7) +
  ylim(0, 7500) +
  xlim(0, 0.01) +
  ggtitle('bootstrap distributions')
ggsave('results/gene-dxy-bootstrap-dens.pdf', width=10, height=12)

TView <- function(d) {
  View(d)
  d
}

# plot smooth and boot for bigchroms
tdxy_boot %>%
  filter(chrom %in% bigchroms) %>%
  select(chrom, zf_pos, boot = boot_q99, smooth) %>%
  gather(type, value, boot:smooth) %>% TView %>%
  ggplot(aes(zf_pos, value, colour = type)) +
  geom_line() +
  facet_wrap(~chrom, ncol = 1)

# try plotting of the ratio
tdxy_boot %>%
  filter(chrom %in% bigchroms) %>%
  select(chrom, zf_pos, boot = boot_q95, smooth) %>%
  ggplot(aes(zf_pos, smooth / boot)) +
  geom_hline(yintercept = 1, alpha = 0.3, colour = "red") +
  geom_line() +
  facet_wrap(~chrom, ncol = 1, strip.position = "left")
ggsave('results/dxy-gene-ratio-95.pdf', width = 8, height = 12)

# smooth / bootstrap ratio with several quantiles
# and a giude at 1 (above means island)
tdxy_boot %>%
  filter(chrom %in% bigchroms) %>%
  select(chrom, zf_pos, smooth, starts_with("boot")) %>%
  gather(quant, boot, starts_with("boot")) %>%
  ggplot(aes(zf_pos, smooth / boot, colour = quant)) +
  geom_hline(yintercept = 1,
             alpha = 0.3,
             colour = "red") +
  geom_line() +
  facet_wrap(~chrom, ncol = 1, strip.position = "left")
ggsave('results/dxy-gene-ratio-quants.pdf',
       width = 8,
       height = 12)

# smooth and bootstraps (several quantiles)
# plotted as points
tdxy_boot %>%
  filter(chrom %in% bigchroms) %>%
  select(chrom, zf_pos, smooth, starts_with("boot")) %>%
  gather(quant, boot, starts_with("boot")) %>%
  ggplot(aes(zf_pos, boot, colour = quant)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_line(aes(y = smooth, group = chrom),
            data = tdxy_boot %>% filter(chrom %in% bigchroms),
            size = 2, colour = "gray70",
            alpha = 0.7) +
  facet_wrap(~chrom, ncol = 1, strip.position = "left") +
  ylim(0, 0.010)
ggsave('results/gene-dxy-bootstraps-point.pdf',
       width = 10,
       height = 16)

#### convert to islands ####
source("interval-tools.R")

tdxy_boot %>%
  filter(smooth > boot_max) %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_max

tdxy_boot %>%
  filter(smooth > boot_q99) %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_99

tdxy_boot %>%
  filter(smooth > boot_q75) %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_75

tdxy_boot %>%
  filter(smooth > boot_q95) %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_95

#### genes in islands ####

read_tsv('results/oo-genes-mart-wiki.txt') ->
  oo_genes

# chceck the overlap of oo genes and windows manually
# - this does not work well, genes are to small;)
oo_genes %>%
  select(chrom, start=start_position, end=end_position) %>%
  mutate(type="oo genes",
         chrom=as.character(chrom)) %>%
  bind_rows(d_wins_99 %>% mutate(type="win 99")) %>%
  mutate(id = row_number()) %>%
  gather(ptype, pos, start, end) %>%
  ggplot(aes(pos, y=type, , colour=type, group=id)) +
  geom_line(size=2) +
  facet_wrap(~chrom, ncol = 1, switch = "y")

d_wins_99 %>%
  mutate(id = row_number(),
         type = "windows") %>%
  gather(ptype, pos, start, end) %>%
  ggplot(aes(pos, y = type, colour = type)) +
  geom_line(aes(group = id), size = 2) +
  geom_point(aes(x = start_position),
             data = oo_genes %>%
               mutate(chrom = paste0("chr", chrom), type = "genes")) +
  facet_wrap( ~ chrom, ncol = 1, strip.position = "left")
ggsave('results/gene-dxy-oo-genes-q95.pdf',
       width = 10,
       height = 16)


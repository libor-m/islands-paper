#### song related probes ####
#
# name the importatn columns
colnames(blat)[c(10, 14, 16, 17)] <- c("probe", "chrom", "zf_start", "zf_end")

# annotate by the original module annotation
probes <- read.delim("song/song-modules.tsv", header=F) %>%
  dplyr::rename(module=V6) %>%
  dplyr::select(module) %>%
  mutate(probe=paste0("seq", 1:length(module)))

blat <- read.delim('song/probe-blat-full.psl', skip=5, header=F) %>%
  dplyr::rename(match=V1, mismatch=V2, probe=V10, chrom=V14, zf_start=V16, zf_end=V17) %>%
  filter(match > 58, mismatch < 2) %>%
  mutate(zf_pos = (zf_start + zf_end) / 2) %>%
  inner_join(probes)

irblat <- blat %>% {IRanges(start=.$zf_start, end=.$zf_end)}
grblat <- blat %>% {GRanges(seqnames=.$chrom, ranges=irblat)}
igblat <- findOverlaps(gra, grblat) %>%
  as.data.frame %>%
  mutate(chrom = as.factor(seqnames(grblat)[subjectHits]),
         zf_pos = blat %>% {.$zf_pos[subjectHits]},
         module = blat %>% {.$module[subjectHits]})

# TODO: should split the plot into premade parts not to copy the code around
tm %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) +
  geom_line(aes(y=fst_boot), colour="yellow") +
  geom_line(aes(y=fst_smooth), colour="blue") +
  geom_point(y=.2, colour="blue", shape=15, size=2, data=tm %>% filter(fst_smooth > fst_boot, chrom %in% bigchroms)) +
  geom_point(aes(shape=module), data=blat %>% filter(chrom %in% bigchroms), size=3, colour="red", y=-.1) +
  scale_shape_manual(values=c(1, 2, 3)) +
  facet_wrap(~chrom, ncol=1) +
  ylim(c(-.1, .2)) +
  ggtitle("nightingale speciation islands as mapped to zebra finch chromosomes, 25k bootstrap, song gene probes")

ggsave('results/fst_islands-song.pdf', width=20, height=16)
ggsave('results/fst_islands-song-all-cand.pdf', width=20, height=16)

ggplot(blat, aes(V1)) + geom_histogram()


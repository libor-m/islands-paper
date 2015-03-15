library(dplyr)
library(ggplot)
library(gtools)

# read the data and remove nan values
d <- read.delim("data/lp2-var-filtered.weir.fst") %>% filter(!is.na(WEIR_AND_COCKERHAM_FST))

# check the maximum sizes of chromosomes
clens <- d %>% 
  group_by(CHROM) %>%
  summarise(len=max(POS)) %>%
  arrange(desc(len))

# rearrange the factor and plot chromosome lengths
clens %>%
  mutate(CHROM=factor(CHROM, levels=rev(CHROM))) %>%
  ggplot(aes(CHROM, len)) + geom_bar(stat="identity") + coord_flip()

# pick first n biggest chromosomes
n <- 10
d %>%
  filter(CHROM %in% clens$CHROM[1:n]) %>%
  mutate(CHROM=factor(CHROM, levels=CHROM %>% levels %>% mixedsort)) %>%
  ggplot(aes(POS, WEIR_AND_COCKERHAM_FST)) + 
    geom_point() +
    facet_wrap(~CHROM, nrow=n)

# calculate rolling mean Fst
# equals to 'box' kernel
win <- 10000
dc <- d %>%
  mutate(w=cut(POS, seq(0, max(POS)+win+1, by=win)))
  group_by(CHROM, w) %>%
  summarise(POS=mean(POS), WEIR_AND_COCKERHAM_FST=)
library(tidyverse)

# Variants and Fst values ----

# change chrom column in df to a factor with levels
# correctly sorted for plotting with ggplot

dvar <- read.delim("data/variant-table.tsv")
dfst <- read.delim('data-genome/lp2-var-filtered.weir.fst', col.names=c("chrom", "pos", "fst"))

# merge fst into variants information
# some rows in dfst are missing the value, filter those out
# filter out contigs badly mapped to zebra finch reference
# filter out low quality variants
dvar %>%
  inner_join(dfst) %>%
  filter(!is.na(fst)) %>%
  filter(zf_max - zf_min < 5e5) %>%
  filter(qual > 10) ->
  daf

# add more specific name for interop with dxy
daf_fst <- daf

# save filtered data
save(daf_fst, file='data/daf_fst.RData')
write.table(daf, file='data/daf-fst.tsv', row.names=F, quote=F, sep="\t")

# load when resuming session
load('data/daf_fst.RData')

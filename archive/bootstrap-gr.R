#
# functions to do bootstrapping for sparse float values
# (chromosome "chr", position "zf_pos", value in argument controllable column
#

#
# the idea is to find all variants for each (overlapping)
# window once, store the index to the variant table
#
# the bootstrap then shuffles the [Fst] values in variant table
# and calculates the smoothed value (aggregated by max)
#
# the real data is obtained with the same smoothing, without shuffling
#
# this shadows half of R base functions and clashes with dplyr (wtf..)
# TODO: consider using data.table::foverlaps instead of this behemoth
library(GenomicRanges)

# this has to be done only once, subsequent queries
# with bootstrapped fst can use the result
# vars is a df with chrom, zf_pos, zf_max, zf_min
find_variants <- function(vars, win_size=1e6, req_points=1e3) {
  # ensure that the whole biggest chrom is covered with windows
  stride <- as.integer(max(vars$zf_max) %>% { .  / max(. / win_size, req_points) })
  
  # create equally spaced poins along all chromosomes, spaced by `stride`
  dfpoints <- vars %>% 
    filter(chrom != "chrUn") %>%         # take only the known chromosomes
    mutate(chrom=as.character(chrom)) %>% # do not copy around the factor levels
    group_by(chrom) %>%                   # calculate chromosome sizes
    summarize(zf_size=max(zf_max)) %>%    # ..
    filter(zf_size > win_size) %>%        # pick chromosomes bigger than requested window
    split(seq_len(dim(.)[1])) %>%         # make data suitable for lapply
    lapply(function(x) 
      data.frame(chrom=x[[1]], zf_pos=seq(from=win_size/2, to=x[[2]], by=stride), stringsAsFactors=F)) %>%
    rbind_all                         # bind the dataframes together
  
  # place variants without direct mapping in zebra finch into wide range of all exons in given contig
  ir <- IRanges(start=ifelse(is.na(vars$zf_pos), vars$zf_min, vars$zf_pos), 
                end=ifelse(is.na(vars$zf_pos), vars$zf_max, vars$zf_pos + 1))
  
  # create interval forest for faster querying
  gr <- GRanges(seqnames=vars$chrom, ranges=ir)
  gf <- GIntervalTree(gr)
  
  # create whole set of query ranges
  # windows are centered on the equally spaced points
  irq <- IRanges(start=dfpoints$zf_pos - win_size / 2, end=dfpoints$zf_pos + win_size / 2)
  gq <- GRanges(seqnames=dfpoints$chrom, ranges=irq)
  
  # query the intervalforest
  # annotate the results with chromosome and midpoint of query window
  rv <- findOverlaps(gq, gf) %>% 
    as.data.frame %>%
    mutate(chrom=as.factor(seqnames(gf)[subjectHits]), 
           zf_pos=gq %>% ranges %>% start %>% {.[queryHits] + win_size / 2}  )
  
  rv
}

smoothed_values <- function(hits, values) {
  hits %>%
    mutate(fst=values[subjectHits]) %>%
    group_by(chrom, zf_pos) %>%
    summarize(smooth=mean(fst, na.rm=T), nvars=n()) %>%
    ungroup
}

# `sample()` is insane, sampling a vector of length 1 - c(N)
# returns a random permutation of ints < N
sane_sample <- function(x) {
  if(length(x) > 1) {
    sample(x)
  } else {
    x
  }
}


# fst randomized by sampling the same chromosome
rand_fst <- function(d)
  d %>% 
  select(chrom, zf_pos, fst) %>%
  group_by(chrom) %>%
  mutate(fst=sane_sample(fst)) %>%
  .[,"fst"] %>%
  .[[1]]

# random permutation in the same chromosome
# mutate instead of sample alone to permute inside chromosomes
# .dots is a bit excersising of NSE, the resulting name could be fixed..

# FUCK R .. spent an hour debugging this:
# If x has length 1, is numeric (in the sense of is.numeric) and x >= 1, 
# sampling via sample takes place from 1:x. 
# Note that this convenience feature may lead to undesired behaviour 
# when x is of varying length in calls such as sample(x). 
# .. sample()'s surprise -- example
rand_var0 <- function(d, var="fst") {
  sane_sample <- function(x) {
    if(length(x) > 1) {
      sample(x)
    } else {
      x
    }
  }
  
  d %>% 
  select_("chrom", "zf_pos", var) %>%
  group_by(chrom) %>%
  mutate_(rand=paste0("sane_sample(", var, ")")) %>%
  .[,"rand"] %>%
  .[[1]]
}

# need dplyr v0.4+
# https://github.com/hadley/dplyr/issues/708
rand_var <- function(d, col="fst") {
  d %>%
    select_("chrom", col) %>%
    rename_(.dots = setNames(col, "var")) %>%
    mutate(rand = sane_sample(var)) %>%
    .[,"rand"] %>%
    .[[1]]
}

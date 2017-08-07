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
library(tidyverse)

library(data.table)

find_variants <- function(vars, window_size=1e6, stride=1e5) {

    # calculate 'window' coordinates along all chromosomes
  vars %>% 
    filter(chrom != "chrUn") %>%             # take only the known chromosomes
    mutate(chrom = as.character(chrom)) %>%  # do not copy around the factor levels
    group_by(chrom) %>%                      # calculate chromosome sizes
    summarize(                               # find range of var coords for each chrom
      zf_min = min(zf_min),
      zf_max = max(zf_max)) %>% 
    filter(zf_max > window_size) %>%         # pick chromosomes bigger than requested window
    purrrlyr::by_row(
      .collate = "rows",
      .to = "start",
      function(x)
      seq(from = 0,
          to = x$zf_max,
          by = stride)) %>% 
    select(chrom, start) %>% 
    mutate(end = start + window_size) %>%
    setDT %>%
    setkey(chrom, start, end) ->
    dt_wins

  # insert variants with unsure mapping to genome
  # into interval of surrounding exons
  # (ie gene is mapped, but the particular exon/intron with the variant 
  # is not)
  vars %>%
    mutate(start = coalesce(zf_pos, zf_min),
           end = coalesce(zf_pos + 1, zf_max),
           src_row = row_number()) %>%
    select(chrom, start, end, src_row) %>% 
    setDT ->
    dt_vars
    
  # for each interval find variants which fall into it
  # dt_wins is small table with large spans
  # dt_vars is big table
  foverlaps(dt_vars, dt_wins) %>%
    mutate(window_midpoint = start + window_size / 2) -> 
    dt_hits
}

smoothed_values <- function(hits, values) {
  hits %>%
    mutate(fst = values[src_row]) %>%
    group_by(chrom, window_midpoint) %>%
    summarize(smooth = mean(fst, na.rm = T), 
              nvars = n()) %>%
    rename(zf_pos = window_midpoint) %>%   # legacy compatibility ;(
    ungroup
}

# `sample()` is insane, sampling a vector of length 1 
# `samle(10) == sample(c(10))` random permutation of ints 1:10
sane_sample <- function(x) {
  if (length(x) > 1) {
    sample(x)
  } else {
    x
  }
}

# randomize, using the grouping set outside
rand_var <- function(d, col="fst") {
  d %>%
    select_(.dots = c(d %>% groups, col)) %>%
    rename_(.dots = setNames(col, "var")) %>%
    mutate(rand = sane_sample(var)) %>%
    .$rand
}

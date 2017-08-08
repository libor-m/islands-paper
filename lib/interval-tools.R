#
# a substitution for GRanges, which upon loading
# breaks most of my normal R..
#
#
# try to calculate 'reduce' without granges
#
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

#### implementation #####

# reduce without using GenomicRanges
# margin=1 connects book-ended intervals,
# increasing it merges intervals further apart
reduce_intervals_zf <- function(dints, margin=1)
  dints %>%
    group_by(chrom) %>%
    arrange(zf_start) %>%
    mutate(prev_end = zf_end %>% lag,
           tail = prev_end %>% replace(is.na(prev_end), 0) %>% cummax,
           sw = zf_start > tail + margin,
           newfeat = sw %>% cumsum %>% as.factor) %>% 
    group_by(chrom, newfeat) %>%
    summarise(start = zf_start %>% min,
              end = zf_end %>% max) %>%
    select(-newfeat) %>%
    ungroup

# convert points to windows of win_size centered around each point
points_to_wins <- function(dpoints, win_size=1e6, margin=1)
  dpoints %>% 
    mutate(zf_start=pmax(zf_pos - win_size / 2, 0), 
           zf_end=zf_pos + win_size / 2) %>%
    select(chrom, zf_start, zf_end) %>%
    filter(!is.na(zf_start)) %>%
    reduce_intervals_zf(margin)

# overlap intersection
# found in 
# http://stackoverflow.com/questions/27574775/..
# ../is-it-possible-to-use-the-r-data-table-funcion-foverlaps-to-find-the-intersectio
dt.intersect <- function(query_bigger, subj_smaller) {
  setkey(setDT(query_bigger), chrom, start, end)
  setkey(setDT(subj_smaller), chrom, start, end)

  resultTable <- foverlaps(query_bigger, subj_smaller, nomatch = 0)

  resultTable[, start := pmax(start, i.start)]
  resultTable[, end := pmin(end, i.end)]
  resultTable[, `:=`(i.start = NULL, i.end = NULL)]
  
  as.data.frame(resultTable)
}

# input is chrom, start, end
reduce_intervals <- function(dints, margin=1)
  dints %>%
  group_by(chrom) %>%
  arrange(start) %>%
  mutate(prev_end = end %>% lag,
         tail = prev_end %>% replace(is.na(prev_end), 0) %>% cummax,
         sw = start > tail + margin,
         newfeat = sw %>% cumsum %>% as.factor) %>% 
  group_by(chrom, newfeat) %>%
  summarise(start = start %>% min,
            end = end %>% max) %>%
  select(-newfeat) %>%
  ungroup

# intervals are chrom, start, end
# genome is chrom, start, end
# 1 based coordinates
interval.complement <- function(intervals, genome) {
  intervals %>%
    reduce_intervals %>% 
    left_join(genome %>% select(chrom, len), by="chrom") %>%
    group_by(chrom) %>%
    do(data.frame(
         start=c(1, .$end),
         end=c(.$start, .$len %>% head(1))
       )) %>% 
    reduce_intervals %>%
    filter(start <= end)
}

#### dev code ####
if(FALSE) {
win_size <- 1e6
tdxy_boot %>% 
  filter(smooth > boot) %>%
  mutate(zf_start=pmax(zf_pos - win_size / 2, 0), zf_end=zf_pos + win_size / 2) %>%
  select(chrom, zf_start, zf_end) ->
  islands

# use tail + 1 to merge book-ended features
islands %>%
  group_by(chrom) %>%
  arrange(zf_start) %>%
  mutate(prev_end = zf_end %>% lag,
         tail = prev_end %>% replace(is.na(prev_end), 0) %>% cummax,
         sw = zf_start > tail + 1,
         newfeat = sw %>% cumsum %>% as.factor) %>% 
  group_by(chrom, newfeat) %>%
  summarise(zf_start = zf_start %>% min,
            zf_end = zf_end %>% max) %>%
  select(-newfeat) ->
  reduced

# check reduction with plot
reduced %>%
  mutate(id=row_number()) %>%
  gather(type, pos, zf_start, zf_end) ->
  reduced_long

islands %>%
  group_by(chrom) %>%
  mutate(id=row_number()) %>%
  gather(type, pos, zf_start, zf_end) %>%
  ggplot(aes(pos, group=id)) +
  geom_line(aes(y=id)) +
  geom_line(data=reduced_long, y=1.2, size=2) +
  facet_wrap(~chrom, scale="free")
}
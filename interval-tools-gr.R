#
# operation on intervals using genomic ranges 
#

# convert GenomicRanges object to data frame
gr_to_bed3 <- function(gr) {
  data.frame(chrom=gr %>% seqnames, start=gr %>% start, end=gr %>% end)
}

# data frame (with zf coords) to GRanges
points_to_gr <- function(d, win_size) {
  d %>% 
    mutate(zf_start=pmax(zf_pos - win_size / 2, 0), zf_end=zf_pos + win_size / 2) %>%
    select(chrom, zf_start, zf_end) ->
    islands
  
  islands %>% 
{IRanges(start=.$zf_start, end=.$zf_end)} ->
  ir

GRanges(seqnames=islands$chrom, ranges=ir)
}

# convert selected points along chromosomes
# into windows, resolve the overlaps
points_to_wins <- function(d, win_size=1e6) {
  
  d %>%
    points_to_gr(win_size) %>%
    reduce ->
    grr
  
  gr_to_bed3(grr)
}

# convert data frame to GRanes object
win_to_gr <- function(d) {
  GRanges(seqnames=d$chrom, ranges=IRanges(start=d$start, end=d$end))
}

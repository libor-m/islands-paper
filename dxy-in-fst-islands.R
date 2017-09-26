# 
# mapping per contig dxy into fst islands ----
#

library(tidyverse)
# library(data.table)

source('lib/gene-set-enrichment.R')

# read fst windows
read_tsv('data/d_wins.tsv') %>%
  filter(measure == "Fst") -> d_wins

# paper revision 0
read_tsv('data-gene-dxy/contigs_dxy_avg_to_fst_isl.txt') %>%
  rename(chrom = zf_contig_chrom,
         start = zf_contig_start,
         end = zf_contig_end) %>%
  mutate(dxy_rel = dxy_abs / contig_length) ->
  ddxy

# paper revision 1
read_tsv('data-gene-dxy/dxy_abs_by_contig.tab') %>%
  rename(
    chrom = ZF_CONTIG_CHROM,
    start = ZF_CONTIG_START,
    end = ZF_CONTIG_END,
    contig_id = CONTIG,
    contig_length = SIZE,
    dxy_abs = DXY_ABS,
    no_overlapping_islands = FST_ISLAND) %>%
  mutate(dxy_rel = dxy_abs / contig_length) ->
  ddxy

# for each contig find an island that contains the contig in d_wins
# and add it's id to the contig

# d_wins is small table with large spans (-> y, keyed)
# ddxy is big table -> x
# reach into different universe (of data.table) for the awesome
# `foverlaps` function
data.table::setkey(data.table::setDT(d_wins), chrom, start, end)

data.table::foverlaps(data.table::setDT(ddxy), d_wins, type = "within") ->
  d_join

d_join %>%
  filter(!is.na(measure)) -> 
  dovr

# plot intra vs inter fst island dxy
d_join %>%
  mutate(`In Fst Island` = !is.na(measure)) %>%
  ggplot(aes(dxy_rel, fill = `In Fst Island`)) +
  geom_density(alpha = 0.5, colour=NA) +
  xlim(0, 0.01)
ggsave('results/dxy-gene-in-out-islands.pdf', width=8, height=6)

d_join %>%
  ggplot(aes(is.na(measure), dxy_rel)) +
  geom_boxplot()

d_join %>%
  group_by(island = !is.na(measure)) %>%
  summarise(mean = mean(dxy_rel), median = median(dxy_rel))

d_join %>%
  mutate(island = !is.na(measure)) %>%
  select(island, dxy_rel) %>%
  t.test(dxy_rel ~ island, .)

dovr %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  dim

dovr %>% write.table(
  "data-gene-dxy/contigs-in-fst-wins.tsv",
  sep = "\t",
  quote = F,
  row.names = F
)
read_tsv("data-gene-dxy/contigs-in-fst-wins.tsv") -> dovr

# mean relative  dxy per island
dovr %>%
  group_by(chrom, start, end, id) %>%
  summarise(dxy_rel_mean = dxy_rel %>% mean) %>%
  ungroup ->
  dovrmean

dovrmean %>%
  ggplot(aes(dxy_rel_mean)) +
  geom_histogram()

dovrmean %>%
  arrange(desc(dxy_rel_mean)) %>% 
  View

# find some limit for picking the islands
dovrmean %>%  
  .$dxy_rel_mean %>%
  quantile(probs = c(.5, .75, .95))

# mean looks usable as a limit
dovrmean %>%
  .$dxy_rel_mean %>%
  mean ->
  mean_island_dxy

d_join %>%
  .$dxy_rel %>%
  mean ->
  mean_genome_dxy

# get the genes in selected windows ----
source("lib/mart.R")

dovrmean %>%
  filter(dxy_rel_mean > mean_genome_dxy) %>%
  encode_wins %>%
  genes_from_regions ->
  int_genes

# save data for vasek
dovrmean %>%
  mutate(high_dxy = dxy_rel_mean > mean_genome_dxy) %>%
  write.table("data-gene-dxy/fst-islands-mean-dxy.tsv", sep="\t", row.names=F, quote=F)

int_genes %>%
  write.table("data-gene-dxy/fst-high-dxy-islands-ensgenes.tsv", sep="\t", row.names=F, quote=F)

# do the tests ----

source('lib/gene-set-enrichment.R')

# kegg <- kegg_get_data('tgu', 'ensembl-tgu')
kegg <- kegg_get_data_offline()
read_tsv('data-annot/uni_go.tsv') -> uni_go

# kegg$pathway_genid %>% write.table("data-annot/pathway_ensgene.tsv", quote=F, row.names=F, sep="\t")

# here the question is what to use as a background
# - the whole universe: uni_go: uni_go$ensembl_gene_id
#   this gives some interesting results, but without statistical
#   significance after FDR correction ..
# - all genes annotated with a pathway: kegg$pathway_genid$genid
#   this way all the groups seem severely underrepresented
test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  View

test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  filter(p.fisher <= 0.05) %>%
  as.data.frame %>%
  format(digits = 3) %>%
  write.table("results/islands-fst99-dxy-mean-kegg-rev1.txt", sep = "\t", row.names = F, quote = F)

# test enrichment for fst islands only
d_wins %>%
  encode_wins %>%
  genes_from_regions ->
  int_genes

test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  View

# directly save the results
test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  filter(p.fisher <= 0.05) %>%
  format(digits=3) %>%
  write.table("results/KEGG-fst-islands-1M.txt", sep="\t", row.names=F, quote=F)

# more data for vasek
int_genes %>%
  write.table("data-gene-dxy/fst-islands-ensgenes.tsv", sep="\t", row.names=F, quote=F)

read_tsv('data-gene-dxy/fst-high-dxy-islands-ensgenes.tsv') -> int_high
setdiff(int_genes, int_high) -> int_low


# genes in Fst islands with lower Dxy ----

dovrmean %>%
  filter(dxy_rel_mean < mean_genome_dxy) %>%
  encode_wins %>%
  genes_from_regions ->
  int_genes

int_genes %>%
  write.table("data-gene-dxy/fst-low-dxy-islands-ensgenes.tsv", sep="\t", row.names=F, quote=F)

test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  View

# directly save the results
test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  filter(p.fisher <= 0.05) %>%
  as.data.frame %>%
  format(digits = 3) %>%
  write.table("results/KEGG-fst99-dxy-lt-mean-1M.txt", 
              sep = "\t", 
              row.names = F, 
              quote = F)

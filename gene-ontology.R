#
# analyze GO/pathways in detected islands
#

library(tidyverse)
library(gtools)

#### load data ####

tfst_boot <- read_tsv("data/tfst_boot.tsv")
tdxy_boot <- read_tsv("data/tdxy_boot.tsv")

t_boot <- bind_rows(tfst_boot, tdxy_boot)

#### define islands ####

#
# island is defined as an 1M window around any variant, that
# exceeds the bootstrap threshold (the window size is the same used to
# calculate the smooth value)
#

source("lib/interval-tools.R")

# need to convert each measure separately
# because the windows are merged together and there is no 'measure' grouping
tfst_boot %>%
  filter(smooth > boot_q99,
         measure == "Fst") %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_fst

tdxy_boot %>%
  filter(smooth > boot_q95,
         measure == "Dxy") %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_dxy

# intersect windows of elevation for both measures
d_wins_int <- dt.intersect(d_wins_dxy, d_wins_fst)

bind_rows(d_wins_dxy %>% mutate(measure = "Dxy"),
          d_wins_fst %>% mutate(measure = "Fst"),
          d_wins_int %>% mutate(measure = "both")) %>%
  mutate(id = row_number()) ->
  d_wins

# d_wins %>% write.table("data/d_wins.tsv", quote = F, row.names=F, sep="\t")
d_wins %>% write_tsv("data/d_wins.tsv")
read_tsv('data/d_wins.tsv') -> d_wins

# visual check if the intersection is ok
d_wins %>%
  filter(chrom %in% bigchroms) %>%
  gather(type, pos, start, end) %>%
  ggplot(aes(pos, colour=measure)) +
  geom_line(aes(y=measure, group=id), size=3) +
  facet_wrap(~chrom, ncol = 1, strip.position = "left")
ggsave("results/window-overlaps.pdf", width=297, height=210, units="mm")


#### get genes in islands ####
#
# pull chicken homologs from biomart for those dbs, that don't know zebra finch
# all GO terms in finch are IEA anyways
# (http://geneontology.org/page/automatically-assigned-evidence-codes)
source("lib/mart.R")

# convert wins to regions string
# and pull the data
d_wins %>%
  filter(measure == "both") %>%
  encode_wins %>%
  genes_from_regions ->
  int_genes

# test only Fst windows
d_wins %>%
  filter(measure == "Fst") %>%
  encode_wins %>%
  genes_from_regions ->
  int_genes

# check what we got
table(int_genes$ggallus_homolog_orthology_type)
# ~1k one to one orthologs in chicken out of 1377 hits

int_genes %>% write_tsv("data/islands-intersect-r1.genes.tsv")


#### get zebra finch gene universe ####

# download zebra finch 'gene universe'
uni_go <- getBM(c("ensembl_gene_id",
                  "ensembl_transcript_id",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "go_id",
                  "go_linkage_type"), "", "", mart)

write.table(uni_go, file='data/uni_go.tsv', row.names=F, sep="\t", quote=F)
read_tsv('data-annot/uni_go.tsv') -> uni_go


#### custom GO enrichment test ####

source('lib/gene-set-enrichment.R')

# kegg <- kegg_get_data('tgu', 'ensembl-tgu')
kegg <- kegg_get_data_offline()

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
  filter(p.fisher <= 0.05) ->
  kegg_enrich_universe

kegg_enrich_universe %>%
  as.data.frame %>%
  format(digits = 3) %>%
  write.table(
    "results/islands-fst99-dxy95-kegg-rev1.txt",
    sep = "\t",
    row.names = F,
    quote = F
  )

test_enrichment(unique(int_genes$ensembl_gene_id),
                unique(uni_go$ensembl_gene_id),
                kegg$pathway_genid,
                kegg$pathway_desc) %>%
  arrange(p.fisher) %>%
  filter(p.fisher <= 0.05) %>%
  as.data.frame %>%
  format(digits = 3) %>%
  write.table(
    "results/islands-fst99-kegg-rev1.txt",
    sep = "\t",
    row.names = F,
    quote = F
  )

#### get details on the gene list

# extract the genes from the serialized form
# into one ene per row
# do() is a bit difficult with column naming,
# have to use tiny trick with `colanmes<-`
kegg_enrich_universe %>%
  select(description, qin.list) %>%
  filter(grepl("oocyte", description, ignore.case=T)) %>%
  group_by(description) %>%
  do(data.frame(ensg=strsplit(.$qin.list, ","))) %>%
  `colnames<-`(c("description", "ensg")) ->
  oo_genes

# take the unique genes, and pull some more information on
# them from biomart
getBM(c("ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "go_id",
        "go_linkage_type",
        "name_1006"),
      filters = c("ensembl_gene_id"),
      values = oo_genes$ensg %>% unique,
      mart = mart) ->
  oo_mart

oo_mart %>% write.table(
  "results/oo-genes-mart.txt",
  sep = "\t",
  row.names = F,
  quote = F
)

# leave out GO, use just wikigene
getBM(c("ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "hgnc_symbol",
        "wikigene_name",
        "wikigene_description"),
      filters = c("ensembl_gene_id"),
      values = oo_genes$ensg %>% unique,
      mart = mart) ->
  oo_mart

# output table for the paper
oo_mart %>%
  rename(chrom = chromosome_name) %>%
  sortchrom %>%
  arrange(chrom, start_position) %>%
  write.table(
    "results/oo-genes-mart-wiki.txt",
    sep = "\t",
    row.names = F,
    quote = F
  )

# check chromosome co-localization
oo_mart %>%
  group_by(chromosome_name) %>%
  summarise(ngenes = ensembl_gene_id %>% unique %>% length) %>%
  arrange(desc(ngenes))

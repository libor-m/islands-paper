#
# analyze gene ontology in detected islands
#

library(GenomicRanges)

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)

#### load data ####
setwd("c:/work/lab/slavici-clanek/")

tfst_boot <- read.delim("data/tfst_boot.tsv")
tdxy_boot <- read.delim("data/tdxy_boot.tsv")

t_boot <- bind_rows(tfst_boot, tdxy_boot)

#### define islands ####

#
# island is defined as an 1M window around any variant, that 
# exceeds the bootstrap threshold (the window size is the same used to 
# calculate the smooth value)
#

# old version, for reference
source("interval-tools-gr.R")

# the most interesting should be areas on the intersection of islands
d_wins_int <- BiocGenerics::intersect(win_to_gr(d_wins_dxy), win_to_gr(d_wins_fst)) %>% gr_to_bed3


# non-GRanges version
source("interval-tools.R")

# need to convert each measure separately
# because the windows are merged together and there is no 'measure' grouping
t_boot %>%
  filter(smooth > boot,
         measure == "Fst") %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_fst

t_boot %>%
  filter(smooth > boot,
         measure == "Dxy") %>%
  points_to_wins(1e6, 1e6) ->
  d_wins_dxy

# use the new version
d_wins_int <- dt.intersect(d_wins_dxy, d_wins_fst)
# results of old and new are the same except for
# BioC clipping the intervals at 1 instead of 0

bind_rows(d_wins_dxy %>% mutate(measure="Dxy"),
          d_wins_fst %>% mutate(measure="Fst"),
          d_wins_int %>% mutate(measure="both")) %>%
  mutate(id=row_number()) ->
  d_wins

d_wins %>% write.table("data/d_wins.tsv", quote=F, row.names=F, sep="\t")

# visual check if the intersection is ok
d_wins %>% 
  filter(chrom %in% bigchroms) %>%
  gather(type, pos, start, end) %>%
  ggplot(aes(pos, colour=measure)) +
  geom_line(aes(y=measure, group=id), size=3) +
  facet_wrap(~chrom, ncol = 1)
ggsave("results/window-overlaps.pdf", width=297, height=210, units="mm")

# the 2M window works a bit better for the intersection

# will rather use external tool for the intersection with all zf annotations
# get the reduced bed3 list 
d_wins_fst %>% write.table(file='data/islands-fst.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")
d_wins_dxy %>% write.table(file='data/islands-dxy.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")
d_wins_int %>% write.table(file='data/islands-int.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")

#### connect to biomart ####

# connect to ensembl's biomart
library(biomaRt)
mart <- useMart("ensembl")

# check what they have
dmart <- listDatasets(mart)
dmart %>% filter(grepl("Tae", description))

# connect to zebra finch
mart <- useMart("ensembl", "tguttata_gene_ensembl")
# check the filter names
listFilters(mart) %>% View

# check attribute names
listAttributes(mart) %>% filter(grepl("go", name, ignore.case=T)) %>% View
listAttributes(mart) %>% filter(grepl("gal", name, ignore.case=T)) %>% View
listAttributes(mart) %>% filter(grepl("name", name, ignore.case=T)) %>% View

# test mart
getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "go_id", "go_linkage_type", "name_1006"), 
      filters = c("ensembl_transcript_id"),
      values = c("ENSTGUT00000002388"),
      mart = mart)

#### get genes in islands ####
#
# pull chicken homologs from biomart for those dbs, that don't know zebra finch
# all GO terms in finch are IEA anyways 
# (http://geneontology.org/page/automatically-assigned-evidence-codes)
mart <- useMart("ensembl", "tguttata_gene_ensembl")

# encode islands into biomart's format
encode_wins <- function(dwins)
  dwins %>%
    tidyr::extract(chrom, c("chrom"), regex="chr([[:alnum:]]+)") %>% 
    mutate(region=paste(chrom, start, end, sep=":")) %>%
    .$region %>%
    paste(collapse=",")

genes_from_regions <- function(regions) {
  # need to query only single 'tab' in mart,
  # so no go terms here..
  getBM(c("ensembl_gene_id", 
          "ensembl_transcript_id", 
          "ggallus_homolog_ensembl_gene", 
          "ggallus_homolog_orthology_type", 
          "ggallus_homolog_orthology_confidence"), 
        c("chromosomal_region"), 
        regions, 
        mart)  
}

# convert wins to regions string 
# and pull the data
d_wins %>% 
  filter (measure == "both") %>%
  encode_wins %>% 
  genes_from_regions ->
  int_genes

# check what we got
table(int_genes$ggallus_homolog_orthology_type)
# ~1k one to one orthologs in chicken out of 1377 hits

int_genes %>% write.table("data/islands-int-2M.genes.tsv", row.names=F, quote=F, sep="\t")
int_genes %>% write.table("data/islands-int-1M-gap1M.genes.tsv", row.names=F, quote=F, sep="\t")

# grab the one2one and copy those to the web service
int_genes %>%
  filter(ggallus_homolog_orthology_type == "ortholog_one2one") %>%
  dplyr::select(ggallus_homolog_ensembl_gene) %>%
  unique %>%
  write.table(file="data/islands-int-2M.ggal.ensg", quote=F, col.names=F, row.names=F)
# 955 unique gene ids

#### GO of genes in islands ####

# GO enrichment tools, each has it's own set of problems..
#
# choose more here:
# https://www.bioconductor.org/packages/release/BiocViews.html#___GeneSetEnrichment
# http://nar.oxfordjournals.org/content/37/1/1/T2.expansion.html
# 
# AMIGO rte http://amigo.geneontology.org/rte 
# - just redirects to PANTHER, 
# - the redirect fails with longer lists
# PANTHER - http://pantherdb.org/tools/compareToRefList.jsp
# - does not have zebra finch
# - fails to resolve ~130 chicken genes
# DAVID http://david.abcc.ncifcrf.gov/
# - recognizes zebra finch, but probably does not use IEA?
# - (the top GO categories have like 12 genes each from the whole set)
# g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi
# - does not have KEGG ids

# the original 'hits' were on KEGG pathways
# but there's no easy way to get KEGG onthology terms for the genes..?
# http://www.kegg.jp/kegg/docs/keggapi.html



# PANTHER on ggal, "Panther GO-Slim BP" without Bonferroni
# shows overrepresentation of few classes, 
# RNA processing, sulfur compounds
# and underrepresentation of developmental and 'housekeeping'
# (cell adhesion) stuff
# see panther-islands-int-2M.ggal.txt 

# g:Profiler accepts genomic ranges directly
d_wins_int %>%
  tidyr::extract(chrom, c("chrom"), regex="chr([[:alnum:]]+)") %>% 
  mutate(region=paste(chrom, start, end, sep=":")) %>%
  dplyr::select(region) ->
  int_regions

int_regions %>% write.table("data/islands-int.regions", col.names=F, row.names=F, quote=F)

# and g:Profliler has a convenient R interface
library(gProfileR)
int_regions %>% 
  gprofiler(
    organism="tguttata",
    region_query=T,
    significant=F,
    max_p_value=0.3,
    correction_method="fdr",
    src_filter=c("KEGG", "GO:BP")) ->
  gprof_test

gprof_test %>% arrange(p.value) %>% View

# save it for funneling to google doc (over spreadsheets;):
gprof_test %>% 
  arrange(p.value) %>% 
  dplyr::select(term.name, 
                p.value, 
                significant, 
                term.size, 
                overlap.size, 
                recall, 
                precision) %>%
  write.table("data/gprof-int-2M.txt", row.names=F, sep="\t")


# load the old data from DAVID
david <- read.delim("data/DAVID_ens_thresh1_E0.2.txt")
# the totals seem to be a bit low..


# GOstats is the first shot
# https://bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf
# looks promising, especially part 1.3 (version April 16, 2015)
# okay, KEGGFrame requires KEGG.db, which is deprecated..

#
# giving up, will write the code myself
# 
source('gene-set-enrichment.R')

kegg <- kegg_get_data('tgu', 'ensembl-tgu')

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
  format(digits=3) %>%
  write.table("results/KEGG-islands-1M.txt", sep="\t", row.names=F, quote=F)

#### check the genes from the oocyte groups ####

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

oo_mart %>% write.table("results/oo-genes-mart.txt", sep="\t", row.names=F, quote=F)

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
  rename(chrom=chromosome_name) %>% 
  sortchrom %>% 
  arrange(chrom, start_position) %>% 
  write.table("results/oo-genes-mart-wiki.txt", sep="\t", row.names=F, quote=F)

# check chromosome co-localization
oo_mart %>% 
  group_by(chromosome_name) %>%
  summarise(ngenes=ensembl_gene_id %>% unique %>% length) %>%
  arrange(desc(ngenes))
  
#### pull genes from bioMart based on Ensembl id ####

ensg <- grep("^ENS", igenes$name)
ig_ens <- igenes[ensg]

# mine mart
ig_go <- getBM(c("ensembl_gene_id", "ensembl_transcript_id", "go_id", "go_linkage_type", "name_1006"), 
      c("ensembl_transcript_id"), ig_ens$name, mart)

# check how many annotated genes
length(ig_ens$name)
# 1787
ig_go %>% filter(!is.na(go_id)) %>% .$ensembl_gene_id %>% unique %>% length
# 1625
# looks like only few not annotated genes

write.table(ig_go, file='data/ig_go.tsv', row.names=F, sep="\t", quote=F)

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

#### GO topn ####

# do a simple topn list analysis for terms in gene islands
ig_go %>%
  group_by(name_1006) %>%
  summarize(n=n_distinct(ensembl_gene_id)) %>%
  arrange(dplyr::desc(n)) %>%
  View
# no good, need contrast against background

# use http://www.ebi.ac.uk/QuickGO to browse
# check GO.db, interested in GO:0007292 female gamete generation and all descendants
library(GO.db)
goterms = unlist(Term(GOTERM))
whmf = grep("oocyte", goterms)
whmf = grep("female", goterms)
whmf = grep("spermatogenesis", goterms)

goterms[whmf] %>% View

# TODO:
# other GOs to check in GO tree, find appropriate top GO and visualize:
# female genitalia morphogenesis (chrZ hits)
# hits for meiosis, spermatogenesis, male, femala
# 
# GO:0007140 male meiosis
# GO:0007143 female meiotic division
# GO:0030237 female sex determination
# GO:0030238 male sex determination
# GO:0007292 female gamete generation
# GO:0048232 male gamete generation
# GO:0030539 male genitalia development
# GO:0030540 female genitalia development

# GOBPOFFSPRING looks like the correct hash table to use
# get all offspring terms of female gamete generation
GOsubtree <- function(goid) {
  c(goid, GOBPOFFSPRING[goid] %>% as.list %>% unlist)
}

# get the intersection of island genes and GO descendants of female gamete generation
# add the parent term to the filter
fgg <- GOsubtree("GO:0007292")
fgg_pos <- ig_go %>%
  filter(go_id %in% fgg) %>% 
  left_join(ig_ens %>% {data.frame(ensembl_transcript_id=.$name, chrom=seqnames(.), zf_pos=start(.) + (end(.) - start(.)) / 2 )})

# annotate more types of candidate genes
# requires c(goid = group, ..)
# returns df with subtrees generated by giod, where each row has particular group
# (the group could be transferred more effectively than with join, but whatever)
candidate_GOs <- function(clist) {
  lids <- lapply(names(clist), function(x) data.frame(parent=x, go_id=GOsubtree(x)))
  rbind_all(lids) %>%
    left_join(data.frame(parent=names(clist), group=clist)) %>%
    left_join(as.data.frame(GOTERM) %>% .[,2:4]) %>%
    unique
}

# pick and annotate some candidate 'speciation' GO groups
# differentiate M/F prezygotyc / postzygotic groups?
# cand_goids <- candidate_GOs(c("GO:0007292" = "female", "GO:0046660" = "female", "GO:0007283" = "male"))

# the real list
# GO:0007140 male meiosis
# GO:0007143 female meiotic division
# GO:0030237 female sex determination
# GO:0030238 male sex determination
# GO:0007292 female gamete generation
# GO:0048232 male gamete generation
# GO:0030539 male genitalia development
# GO:0030540 female genitalia development
cand_goids <- candidate_GOs(c(
  "GO:0007140" = "male",
  "GO:0007143" = "female",
  "GO:0030237" = "female",
  "GO:0030238" = "male",
  "GO:0007292" = "female",
  "GO:0048232" = "male",
  "GO:0030539" = "male",
  "GO:0030540" = "female"))

# check the overlap of child terms
cand_goids$go_id %>% unique %>% length
# looks all children are unique

# transform GRanges to data frame
# (only chrom, pos, gene name)
gr2df <- function(gr) gr %>% {data.frame(
  ensembl_transcript_id = .$name, 
  chrom = seqnames(.), 
  zf_pos=(start(.) + end(.)) / 2)}

# pick the genes that fall into the candidate GO subtrees
# this is good to check if there is any potential reason for 
# the high fst
ig_cands <- gr2df(ig_ens) %>% 
  inner_join(ig_go %>% filter(go_id %in% bpids)) %>% 
  inner_join(cand_goids)

#TODO: pick candidates from universe instead of islands
# (download also coordinates from mart to uni_go)
# this is good for what? chromosome enrichment maybe
ig_cands <- uni_go %>% 
  mutate(chrom = paste0("chr", chromosome_name), 
         zf_pos = (start_position + end_position) / 2) %>%
  inner_join(cand_goids)

# plot it all
# - fst smoothed and bootstrap
# - putative islands
# - candidate genes for male and female
tm %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=fst_boot), colour="yellow") +
  geom_line(aes(y=fst_smooth), colour="blue") +  
  geom_point(y=.2, colour="blue", shape=15, size=2, data=tm %>% filter(fst_smooth > fst_boot, chrom %in% bigchroms)) +
  geom_point(aes(y=ifelse(group == "male", -.1, -.05), shape=group), data=ig_cands %>% filter(chrom %in% bigchroms), size=3, colour="red") +
  scale_shape_manual(values=c(1, 2)) +
  facet_wrap(~chrom, ncol=1) +
  ylim(c(-.1, .2)) +
  ggtitle("nightingale speciation islands as mapped to zebra finch chromosomes, 25k bootstrap")

ggsave('results/fst_islands.pdf', width=20, height=16)

# now look at GO/BP of the genes in the highest scoring islands
# tm has the information on smoothed ans bootstrapped fst
# find genes for each island and join the fst values
ig_tree <- GIntervalTree(ig_ens)
gra <- islands %>% {GRanges(seqnames=.$chrom, ranges=ir, fst_smooth=.$fst_smooth, fst_boot=.$fst_boot)}
bpids <- Ontology(GOTERM) %>% 
  {data.frame(term=names(.), onto=.)} %>%
  filter(onto == "BP") %>%
  .$term %>%
  as.character

igenes_ann <- findOverlaps(gra, ig_tree) %>%
  as.data.frame %>%
  mutate(chrom = as.factor(seqnames(ig_tree)[subjectHits]), 
         zf_pos = ig_tree %>% start %>% {.[subjectHits]},
         ensembl_transcript_id = ig_tree$name[subjectHits],
         fst_smooth = gra$fst_smooth[queryHits],
         fst_boot = gra$fst_boot[queryHits],
         fst_diff = fst_smooth - fst_boot) %>%
  left_join(ig_go %>% filter(go_id %in% bpids)) %>%
  filter(!is.na(go_id)) %>%
  unique %>%
  arrange(dplyr::desc(fst_diff))   

# save top 1k gene entries in the most differentiated areas
igenes_ann %>% head(1000) %>% write.table(file="data/fst_diff-top-1k.tsv", sep="\t", quote=F, row.names=F)
igenes_ann %>% filter(chrom == "chrZ") %>% write.table(file="data/fst_diff-chrZ.tsv", sep="\t", quote=F, row.names=F)

# could chromosome 4 contain 'song' cluster?
# TODO: find song genes
igenes_ann %>% filter(chrom == "chr4") %>% View

# plot fst difference for chromosomes
islands %>%
  filter(chrom %in% bigchroms) %>%
  mutate(chromo=chrom %>% factor(levels=chrom %>% levels %>% mixedsort)) %>%
  ggplot(aes(chromo, fst_smooth - fst_boot)) + geom_boxplot()

# just to check we do the filtering right
dbpids <- GOTERM %>%
  as.data.frame %>%
  .[,2:4] %>%
  unique %>%
  arrange(go_id) %>%
  filter(Ontology == "BP")


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

#### spare parts ####

# after filtering with bedtools..
library(rtracklayer)
igenes <- import.bed('data/genes-in-islands.bed')

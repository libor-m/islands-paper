#
# calculate gene set enrichment
# input:
#  gene list
#  background list
#  term - gene map (data frame, two columns, multiple rows for single term)
# output:
#  table row for each term, set sizes, intersection sizes, p-values
# 
# method:
#  see the end of doc/Gene_list_enrichment_Mar10.pdf for p-values
#  the rest is just some data wrangling

library(tidyr)
library(dplyr)
#### enrichment test implementation ####
# query and background are vectors
# term_to_gene is a two column multimap (pathway, ensg; repeating same term for all genes)
# term_to_desc is a two column map (pathway, *any*; term description)
test_enrichment <- function(query, background, term_to_gene, term_to_desc) {
  
  # return p.value of fisher's exact test
  # calculate values in the contingency table from the marginals
  # to keep the same interface as vec.binom
  mx.fisher <- function(qin, qlen, bin, blen)
    matrix(data=c(qin, qlen - qin, bin, blen - bin), nrow=2) %>%
    fisher.test %>%
    .$p.value
  
  # return p.value of binomial test
  vec.binom <- function(qin, qlen, bin, blen)
    binom.test(qin, qlen, p=bin/blen)$p.value
  
  # start with the list (df) of term - gene pairs
  term_to_gene %>%
    
    # for each term
    group_by(pathway) %>%
    
    # calculate the set sizes and intersections
    # with query and background
    summarise(bin  = length(intersect(genid, background)), 
              blen = length(background),
              qin  = length(intersect(genid, query)),
              qlen = length(query)) %>%
    
    # for each summary 
    # calculate the tests
    mutate(qexp = qlen * bin / blen,
           fc   = qin / qexp,
           p.fisher = mapply(mx.fisher, qin, qlen, bin, blen),
           p.binom  = mapply(vec.binom, qin, qlen, bin, blen),
           p.fisher.fdr = p.adjust(p.fisher, "fdr"),
           p.binom.fdr  = p.adjust(p.binom, "fdr")) %>%
    
    # attach legible term descriptions
    left_join(term_to_desc)
}

# pull the needed data from KEGG API
# and linkdb api
# see http://www.genome.jp/dbget/prefix.html for linked_db
kegg_get_data <- function(kegg_db, linked_db) {
  # curl http://rest.kegg.jp/link/pathway/tgu > data/kegg-link-pathway-tgu
  read.table(paste0("http://rest.kegg.jp/link/pathway/", kegg_db), 
             header=F, 
             col.names = c("kegg", "pathway")) ->
    d_org_pathway
  
  # try to get ensembl finch gene id to kegg pathway map
  # got tgu_ensembl.list from http://www.genome.jp/linkdb/
  # maps 7832 kegg ids to 7832 genes
  # remove the prefix with tidyr::extract
  read.table(paste("http://rest.genome.jp/link", linked_db, kegg_db, sep="/"),
             header=F, 
             col.names = c("kegg", "kegg_genid", "type")) %>%
    extract(kegg_genid, c("genid"), ".*:([[:alnum:]]+)") %>% 
    select(-type) ->
    d_org_genid
  
  # construct a direct pathway -> ensg mapping
  inner_join(d_org_pathway, d_org_genid) %>%
    select(-kegg) %>%
    arrange(pathway) %>%
    unique ->
    d_pathway_genid
  
  # get pathway descriptions
  read.table(paste0("http://rest.kegg.jp/list/pathway/", kegg_db),
             header = F,
             col.names = c("pathway", "description"), 
             sep="\t") ->
    d_pathway_desc 
  
  list(pathway_genid=d_pathway_genid, pathway_desc=d_pathway_desc)
}

#### testing/development section ####

# do not run this code when being source()'d
if(FALSE) {

# pull the data from KEGG's api
# http://www.kegg.jp/kegg/docs/keggapi.html

# read the pathways
# curl http://rest.kegg.jp/link/pathway/tgu > data/kegg-link-pathway-tgu
# - maps 4249 genes to 165 pathways
read.table("data/kegg-link-pathway-tgu.list", 
           header=F, 
           col.names = c("tgu", "pathway")) ->
  d_tgu_pathway

# try to get ensembl finch gene id to kegg pathway map
# got tgu_ensembl.list from http://www.genome.jp/linkdb/
# maps 7832 kegg ids to 7832 genes
# remove the prefix with tidyr::extract
read.table("data/kegg-link-tgu-ensembl.list", 
           header=F, 
           col.names = c("tgu", "kegg_ensg", "type")) %>%
  extract(kegg_ensg, c("ensg"), ".*:([[:alnum:]]+)") %>% 
  select(-type) ->
  d_tgu_ensg

# construct a direct pathway -> ensg mapping
inner_join(d_tgu_pathway, d_tgu_ensg) %>%
  select(-tgu) %>%
  arrange(pathway) %>%
  unique ->
  d_pathway_ensg

# get pathway descriptions
read.table("http://rest.kegg.jp/list/pathway/tgu",
           header = F,
           col.names = c("pathway", "description"), 
           sep="\t") ->
  d_pathway_desc 

# load the data from previous steps
d_wins <- read.delim("data/d_wins.tsv")
int_genes <- read.delim("data/islands-int-2M.genes.tsv")
uni_go <- read.delim("data/uni_go.tsv")

# pull the unique ids
query <- unique(int_genes$ensembl_gene_id)
background <- unique(uni_go$ensembl_gene_id)
background <- unique(d_tgu_ensg$ensg)
  
# create a matrix, do a fisher test, return p-value
# calculate values in the contingency table from the marginals
mx.fisher <- function(qin, qlen, bin, blen)
  matrix(data=c(qin, qlen - qin, bin, blen - bin), nrow=2) %>%
  fisher.test %>%
  .$p.value

# p.value of binomial test
vec.binom <- function(qin, qlen, bin, blen)
  binom.test(qin, qlen, p=bin/blen)$p.value

d_pathway_ensg %>%
  group_by(pathway) %>%
  summarise(bin=length(intersect(ensg, background)), 
            blen=length(background),
            qin=length(intersect(ensg, query)),
            qlen=length(query)) %>%
  mutate(qexp=qlen * bin / blen,
         fc=qin / qexp,
         p.fisher=mapply(mx.fisher, qin, qlen, bin, blen),
         p.binom=mapply(vec.binom, qin, qlen, bin, blen)) %>%
  mutate(p.fisher.fdr=p.adjust(p.fisher, "fdr"),
         p.binom.fdr=p.adjust(p.binom, "fdr")) %>%
  left_join(d_pathway_desc) %>%
  arrange(p.fisher) %>%
  # filter(grepl("oo", description, ignore.case=T)) %>%
  View

}
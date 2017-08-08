# an unstructured mess

# will rather use external tool for the intersection with all zf annotations
# get the reduced bed3 list
d_wins_fst %>% write.table(file='data/islands-fst.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")
d_wins_dxy %>% write.table(file='data/islands-dxy.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")
d_wins_int %>% write.table(file='data/islands-int.bed', row.names=F, col.names=F, sep="\t", quote=F, eol="\n")


# grab the one2one and copy those to the web service
int_genes %>%
  filter(ggallus_homolog_orthology_type == "ortholog_one2one") %>%
  dplyr::select(ggallus_homolog_ensembl_gene) %>%
  unique %>%
  write.table(file = "data/islands-int-2M.ggal.ensg", quote=F, col.names=F, row.names=F)
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

library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", "tguttata_gene_ensembl", host="www.ensembl.org")

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

# spare parts
if (FALSE) {
  listMarts(host="www.ensembl.org")
  mart <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  
  # check what they have
  dmart <- listDatasets(mart)
  dmart %>% filter(grepl("Tae", description))
  
  # connect to zebra finch
  mart <- useMart("ENSEMBL_MART_ENSEMBL", "tguttata_gene_ensembl", host="www.ensembl.org")
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
  
}
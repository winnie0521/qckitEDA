annot_data <- function(dat_filt){
  ensembl = biomaRt::useMart("ensembl",dataset="celegans_gene_ensembl")
  filters = biomaRt::listFilters(ensembl)
  attrs.to.get = biomaRt::listAttributes(ensembl)
  common_ids = biomaRt::getBM(attributes="ensembl_gene_id", filters='wormbase_gene',value=dat_filt$wide$gene, mart = ensembl)
  annoDat <- biomaRt::getBM(attributes=c("ensembl_gene_id",attrs.to.get[c(15,6:9,20,21,5),1]),filters = 'ensembl_gene_id',values =common_ids, mart = ensembl)
  row.names(annoDat) <- annoDat$ensembl_gene_id
  return(annoDat)
  }

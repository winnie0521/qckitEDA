#' Create heatmap of the raw CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @param prsp the choice of whether to use pearson correlation or spearson correlation
#' @return the heatmap of raw CPM
#' @example raw_cpm_box(dat_filt,prsp)

raw_cpm_heatmap <- function(dat_filt,prsp){

  tmpDat.cpm <- reshape2::dcast(dat_filt$long,gene~Sample,value.var = "cpm")[,-c(1)]
  geneCor.pr <- cor(tmpDat.cpm, use="complete.obs", method=paste(prsp))
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdGy")))(n = 300)
  png(filename="heatmap_raw_cpm.png")
  gplots::heatmap.2(geneCor.pr, trace="none", dendrogram="column",col=my_palette, scale="none",cexRow=0.5,cexCol = 0.5)
  dev.off()
}

#' Create heatmap of the log CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @param prsp the choice of whether to use pearson correlation or spearson correlation
#' @return the heatmap of log CPM
#' @example log_cpm_heatmap(dat_filt)


log_cpm_heatmap <- function(dat_filt,prsp){
  tmpDat.cpm <- reshape2::dcast(dat_filt$long,gene~Sample,value.var = "cpm")[,-c(1)]
  geneCor.log.pr <- cor(log(tmpDat.cpm+0.000001), use="complete.obs", method=paste(prsp))
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdGy")))(n = 300)
  gplots::heatmap.2(geneCor.log.pr, trace="none", dendrogram="column",col=my_palette, scale="none")
}

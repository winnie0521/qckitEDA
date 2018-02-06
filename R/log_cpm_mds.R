#' Create MDSplot of the log CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @param designpath the path to the file that contains the design matrix
#' @return the MDSplot of log CPM
#' @example log_cpm_mds(dat_filt,designpath)

log_cpm_mds <- function(dat_filt,designpath){
  model.matrix <- read.csv(designpath,header=T)
  dat_filt$long$logcpm <- log(dat_filt$long$cpm+0.00001)
  tmpDat.cpm <- reshape2::dcast(data=dat_filt$long,gene~Sample, value.var ="logcpm")
  g1<- ggMDSplot(tmpDat.cpm,modMat =model.matrix, modCol = 2,sampleLoc = "Sample", txtSize = 3)
  return(g1)
}

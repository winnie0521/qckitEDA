#' Create MDSplot of the raw CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @param designpath the path to the file that has the design matrix
#' @return the MDSplot of raw CPM
#' @example raw_cpm_mds(dat_filt)



raw_cpm_mds <- function(dat_filt,designpath){
  model.matrix <- read.csv(designpath,header=T)
  tmpDat.cpm <- reshape2::dcast(data=dat_filt$long,gene~Sample, value.var ="cpm")
  g1<- ggMDSplot(tmpDat.cpm,modMat =model.matrix, modCol = 2,sampleLoc = "Sample", txtSize = 3)
  print(g1)
}

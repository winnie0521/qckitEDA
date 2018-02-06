#' Create boxplot of the log CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @return the boxplot of log CPM count of each sample
#' @example log_cpm_box(dat_filt)



log_cpm_box <- function(dat_filt){
  g1 <- ggplot(dat_filt$long,aes(x=Sample,y=log(cpm+0.0000001), fill=Gtype))+geom_boxplot()
  g1+theme(axis.text.x=element_text(size=8, angle=90))+xlab("")
  print(g1)
}

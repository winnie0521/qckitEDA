#' Create boxplot of the raw CPM value of sampels
#' 
#' @param dat_filt the pre-processed data of samples
#' @return the boxplot of raw CPM
#' @example raw_cpm_box(dat_filt)

raw_cpm_box <- function(dat_filt){
  g1 <- ggplot(dat_filt$long,aes(x=Sample,y=cpm, fill=Gtype))+geom_boxplot()+scale_y_log10()
  g1 <- g1+theme(axis.text.x=element_text(size=8, angle=90))+xlab("")+ggtitle("Boxplot of raw CPM")
  return(g1)
}

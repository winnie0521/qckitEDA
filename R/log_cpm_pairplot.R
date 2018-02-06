log_cpm_pairplot <- function(dat_filt){
  dat_filt$long$logcpm <- log(dat_filt$long$cpm+0.00001)
  tmpDat.cpm.log <- reshape2::dcast(dat_filt$long,gene~Sample,value.var = "logcpm")[,-c(1)]
  names(tmpDat.cpm.log)
  png(filename="pairplot.png")
  psych::pairs.panels(tmpDat.cpm.log)
  dev.off()
}

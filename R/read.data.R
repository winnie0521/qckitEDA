#' read in the data and reshape it according to design matrix
#'
#' @param path the path to the data file
#' @param designpath the path to the file that contains the design matrix
#' @return the processed data of the original data file

read.data <- function(path,designpath){
  dat <- read.csv(path, header=T)
  rows_to_keep= apply(dat[,-c(1)],1,FUN=function(x)
  {
    y <- 0;
    if (all(x > 0 )) y <-1;
    return(y)
  }
  )

  dat_filt <- list()
  dat_filt[["wide"]]<- dat[as.logical(rows_to_keep > 0),]
  dat_filt[["long"]]<- reshape2::melt(dat_filt$wide, id="gene", variable.name="Sample")
  model.matrix <- read.csv(designpath, header=T)
  dat_filt[["long"]]<- merge(dat_filt$long, model.matrix, by="Sample", all=T)
  libSizes <- aggregate(data=dat_filt$long, value~Sample,FUN=sum)
  names(libSizes)[2] <- "libsize"
  dat_filt$long <- merge(dat_filt$long, libSizes, by="Sample", all=T)
  dat_filt$long$cpm <- dat_filt$long$value*1000000/dat_filt$long$libsize
  dat_filt$long$logcpm <- log(dat_filt$long$cpm+0.00001)
  return(dat_filt)
}

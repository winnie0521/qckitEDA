#' create PCA plot out of sample data
#' @param dat_filt the dataset that is created from the read.data function
#'

PCA <- function(dat_filt,designpath){
  geneLevelData1 <- dat_filt$wide
  rownames(geneLevelData1) = geneLevelData1[,1]
  geneLevelData1<-geneLevelData1[,-1]
  newmat <- as.matrix(geneLevelData1)
  model.matrix <- read.csv(designpath, header=T)
  modelannot <- AnnotatedDataFrame(data.frame(model.matrix[,2],row.names=colnames(geneLevelData1)))
  data <- EDASeq::newSeqExpressionSet(newmat,
                              phenoData=modelannot)
  EDASeq::plotPCA(data,col=rep(1:4, each=3),labels=TRUE,main="PCA plot of sample gene")
}

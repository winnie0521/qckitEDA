bias.baoxplot <- function(){
  data(geneLevelData)
  data(yeastGC)
  sub <- intersect(rownames(geneLevelData), names(yeastGC))
  submat <- as.matrix(geneLevelData[sub,])
  data <- newSeqExpressionSet(mat,
                              phenoData=AnnotatedDataFrame(
                                data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),
                                           row.names=colnames(geneLevelData))),
                              featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))
  biasPlot(data,"gc",ylim=c(0,5),log=TRUE)
}

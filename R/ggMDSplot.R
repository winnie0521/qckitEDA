#' Create MDS plot out of the wide format of the data
#' 
#' @param wideDat the wide format of data
#' @param cols the column index
#' @param modMat 
#' @param modcol
#' @param txtSize the size of text in the plot
#' @param sampleLoc the index of the sample
#' 
#' @return MDSplot

ggMDSplot <- function(wideDat, cols=c(1,2), modMat=mod.Matrix, modCol=3, txtSize=1, sampleLoc="newSampId")
{

  d <- dist(t(wideDat[,-c(1)]))

  mds.fit <- cmdscale(d, eig = TRUE, k = 2)

  mds.d <- data.frame(x1 = mds.fit$points[, cols[1]], x2 = mds.fit$points[, cols[2]],
                      samples = colnames(wideDat[,-c(1)]))
  mds.d$treatment <- as.factor(modMat[match(as.character(mds.d$samples),as.character(modMat[,sampleLoc])),modCol])

  g1 <- ggplot(mds.d) + geom_text(aes(x = x1, y = x2, label =samples , colour = treatment), size=txtSize)
  return(g1)
}

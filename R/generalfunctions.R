myBarPlot <- function(dat,
                      xVal,
                      yVal,
                      fillVal="Stage",
                      ylabel,
                      xlabel=NULL,
                      ymax=1,
                      abline=FALSE,
                      ablinePos=NULL,
                      logy=FALSE,
                      axs.txt.y=8,
                      discrete.y=FALSE)
{
    a <- ggplot(dat,aes_string(x=xVal,y=yVal,fill=fillVal))+ylim(0,ymax)
    if(logy) a <- a+scale_y_log10()
    if(discrete.y) a <- a+scale_y_continuous(breaks=pretty_breaks())

    a <- a+geom_bar(stat="identity", position="dodge")
    a <- a+ylab(ylabel)+labs(fill = "Condition")

    if (!is.null(xlabel))
        {
            print(xlabel)
            a <- a+xlab(xlabel)
        }

    a <- a+theme(panel.grid.minor=element_blank(),
                 panel.grid.major=element_blank(),
                 legend.background = element_rect(colour = "black"),
                 axis.ticks.x=element_blank(),
                 axis.text.x = element_text(size=6,angle=90, vjust=.5),
                 axis.text.y = element_text(size=axs.txt.y),
                 axis.title.y=element_text(size=10),
                 axis.title.x=element_text(size=10))

    ## Print a set of horizontal cut off lines
    if (abline)
        {
            for (i in  1: length(ablinePos))
                {
                    val <- ablinePos[i]
                    a <- a+geom_hline(aes_string(yintercept=val))
                }
        }
    print(a)
}
## Multiple plot function

## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.

## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.

multiplot <- function(...,
                      plotlist=NULL,
                      file,
                      cols=1,
                      layout=NULL)
    {
        require(grid)

                                        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)

        numPlots = length(plots)

                                        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout))
            {
                                        # Make the panel
                                        # ncol: Number of columns of plots
                                        # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
            }

        if (numPlots==1)
            {
                print(plots[[1]])

            }
        else
            {
                                        # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                                        # Make each plot, in the correct location
                for (i in 1:numPlots)
                    {
                                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                              layout.pos.col = matchidx$col))
                    }
            }
    }


createAnalysisData<-function(inFile)
{
                                        ##inFile= "cortexGeneCounts_sense.sum.cnts.csv"
                                        ##inFile = "cerebGeneCounts_sense.sum.cnts.csv"
    inputDat<-read.csv(file = inFile, header=T)
    zerosToRm <- zeroPositions(inputDat[,-c(1:3)])
    inputDat <- inputDat[zerosToRm,]
    inputDatLong <-melt(inputDat,id.vars=c("gene","chr","start","end"))
    names(inputDatLong)[5] <- "SampleID"
    print ("Converting data to Long Format")
    print(head(inputDatLong))
    libSizes <- aggregate(value~SampleID,data=inputDatLong,sum)
    names(libSizes)[2]<-"libSize"
    print ("\nCalculating Library Sizes")
    print(head(libSizes))
    inputDatLong<-merge(inputDatLong, libSizes, by="SampleID")
    inputDatLong$cpm <- inputDatLong$value*1000000/inputDatLong$libSize
    results <- list(long=inputDatLong, wide=inputDat)
    return(results)
}

createAnalysisDataFromDb<-function(inDat)
{
    library(reshape2)
    tmpDat.wide <- dcast(inDat,geneId~id,value.var="expression")
    tmpDat.1 <- subset(inDat, select=c("geneId","location"))
    tmpDat.1 <- tmpDat.1[!duplicated(tmpDat.1[,"geneId"]),]
    tmpDat.wide <- merge(tmpDat.1,tmpDat.wide, by="geneId")
    tmpDat.wide.chr  <- colsplit(tmpDat.wide$location, ":" , names = c("chr","position"))
    tmpDat.wide.pos <- colsplit(tmpDat.wide.chr$position, "-", names=c("start","end"))
    tmpDat.wide <- tmpDat.wide[,-c(2)]
    tmpDat.wide <- cbind(tmpDat.wide.chr[,1], tmpDat.wide.pos[,1],tmpDat.wide.pos[,2],tmpDat.wide)
    names(tmpDat.wide)[1:4] <- c("chr","start","end","gene")

    print ("\nOriginal data\n")
    print(head(tmpDat.wide))

    print ("\nDimensions of Original data\n")
    print(dim(tmpDat.wide))
    print("\n\n")
    inputDat <- tmpDat.wide
    rm(tmpDat.wide)
    zerosToRm <- zeroPositions(inputDat[,-c(1:3)])
    inputDat <- inputDat[zerosToRm,]
    inputDatLong <-melt(inputDat,id.vars=c("gene","chr","start","end"))
    names(inputDatLong)[5] <- "SampleID"

    print ("\nConverting data to Long Format\n")
    print(head(inputDatLong))
    libSizes <- aggregate(value~SampleID,data=inputDatLong,sum)
    names(libSizes)[2]<-"libSize"

    print ("\nCalculating Library Sizes\n")
    print(head(libSizes))
    inputDatLong<-merge(inputDatLong, libSizes, by="SampleID")
    inputDatLong$cpm <- inputDatLong$value*1000000/inputDatLong$libSize
    results <- list(long=inputDatLong, wide=inputDat)
    return(results)
}



createAnalysisData16p<-function(inFile)
{
  #inFile= "cortexGeneCounts_sense.sum.cnts.csv"
  #inFile = "cerebGeneCounts_sense.sum.cnts.csv"
  inputDat<-read.csv(inFile, header=T)
  zerosToRm <- zeroPositions(inputDat[,-c(1:3)])
  inputDat <- inputDat[zerosToRm,]
  inputDatLong <-melt(inputDat,id.vars=c("gene","chr","start","end"))
  names(inputDatLong)[5] <- "SampleID"
  print(head(inputDatLong))
  idCols<-colsplit(inputDatLong$SampleID, "_", names=c("Animal_ID","Tissue"))
  inputDatLong$AnimalID <- idCols$Animal_ID
  inputDatLong$Tissue <- idCols$Tissue
  print(head(inputDatLong))
  libSizes <- aggregate(value~SampleID,data=inputDatLong,sum)
  names(libSizes)[2]<-"libSize"
  print(head(libSizes))
  inputDatLong<-merge(inputDatLong, libSizes, by="SampleID")
  inputDatLong$cpm <- inputDatLong$value*1000000/inputDatLong$libSize
  results <- list(long=inputDatLong, wide=inputDat)
  return(results)
}

## Function to remove rows with zeros across all samples
zeroPositions <- function(inDat)
{
  row_sub <- apply(inDat[-c(1)], 1, function(row)
      {
          y <- TRUE
          if(all(row ==0 ))
              {y <- FALSE}
          return(y)
      })
  print (paste("Number of Non zero rows:",sum(row_sub)))
  print (paste("Number of rows with zeros:",dim(inDat)[1]-sum(row_sub)))
  return(row_sub)
}

## Function to return the positions from a wide data where all samples have
## expression values greater than a given input threshold

threshPositions <- function(inDat, threshold)
    {
        row_sub <- apply(inDat[-c(1)], 1, function(row)
            {
                y <- TRUE
                if(any(row < threshold ))
                    {y <- FALSE}
                return(y)
            })
        print (paste("Number of rows passing Threshold:",sum(row_sub)))
        print (paste("Number of rows under Threshold:",dim(inDat)[1]-sum(row_sub)))
        return(list(rows=row_sub, genes=inDat$gene[row_sub]))
  }


## Function to recreate dataset from bedtools counts data given genes to remove

filterGenes <- function(inputDat,genelist)
{
    rmPos <- which(as.character(genelist) %in% as.character(inputDat$gene))
    genesToRm <- genelist[rmPos]
    rmPos2 <- which(as.character(inputDat$gene) %in% as.character(genesToRm))
    inputDat.rm <-inputDat[-c(rmPos2),]
    inputDatLong <- melt(inputDat.rm,id.vars=c("gene","chr","start","end"))
    names(inputDatLong)[5] <- "SampleID"
    print(head(inputDatLong))
    print(head(inputDatLong))
    libSizes <- aggregate(value~SampleID,data=inputDatLong,sum)
    names(libSizes)[2]<-"libSize"
    print(head(libSizes))
    inputDatLong<-merge(inputDatLong, libSizes, by="SampleID")
    inputDatLong$cpm <- inputDatLong$value*1000000/inputDatLong$libSize
    results <- list(long=inputDatLong, wide=inputDat.rm)
    return(results)
}



## Function to recreate dataset given genes to remove

filterGenes16p <- function(inputDat,genelist)
{
    rmPos <- which(as.character(genelist) %in% as.character(inputDat$gene))
    genesToRm <- genelist[rmPos]
    rmPos2 <- which(as.character(inputDat$gene) %in% as.character(genesToRm))
    inputDat.rm <- inputDat[-c(rmPos2),]
    inputDatLong <- melt(inputDat.rm,id.vars=c("gene","chr","start","end"))
    names(inputDatLong)[5] <- "SampleID"
    print(head(inputDatLong))
    idCols <- colsplit(inputDatLong$SampleID, "_", names=c("Animal_ID","Tissue"))
    inputDatLong$AnimalID <- idCols$Animal_ID
    inputDatLong$Tissue <- idCols$Tissue
    print(head(inputDatLong))
    libSizes <- aggregate(value~SampleID,data=inputDatLong,sum)
    names(libSizes)[2]<-"libSize"
    print(head(libSizes))
    inputDatLong <- merge(inputDatLong, libSizes, by="SampleID")
    inputDatLong$cpm <- inputDatLong$value*1000000/inputDatLong$libSize
    results <- list(long=inputDatLong, wide=inputDat.rm)
    return(results)
}


## MDSPlots using ggplot2
ggMDSplot <- function(wideDat, cols=c(1,2), modMat=mod.Matrix, modCol=3, txtSize=1, sampleLoc="newSampId")
{
    ## wideDat <- geneDat.wide.rmZ
    ## modMat <- mod.Matrix
    ##cols <- c(1,2)
    ## modCol <- 4
    ##wideDat <- tmpDat.cpm
    ##modMat <- mod.matrix.cere
    ##cols <- c(1,2)
    ##modCol <- 7
    ##sampleLoc <- "SampleID"
    d <- dist(t(wideDat[,-c(1)]))

    mds.fit <- cmdscale(d, eig = TRUE, k = 2)

    mds.d <- data.frame(x1 = mds.fit$points[, cols[1]], x2 = mds.fit$points[, cols[2]],
                        samples = colnames(wideDat[,-c(1)]))
    mds.d$treatment <- as.factor(modMat[match(as.character(mds.d$samples),as.character(modMat[,sampleLoc])),modCol])

    g1 <- ggplot(mds.d) + geom_text(aes(x = x1, y = x2, label =samples , colour = treatment), size=txtSize)
    return(g1)
}

## qqplots using ggplot2
qqGGplot <- function(pvector, colsToUse =NULL)
{
    ##pvector <- lmContr$contr.KV.WV$results$pval
    o <- -log10(sort(pvector,decreasing=F))
    e <- -log10( 1:length(o)/length(o) )
    maxAxis <- -log10(min(pvector[which(pvector > 0)],pvector[which(pvector > 0)]))+1.5
    pdat <-data.frame(list(obs=o,theo=e))
    p1 <- ggplot(pdat,aes(x=theo,y=obs))
    if(is.null(colsToUse)){ p1 <- p1+geom_point()+xlim(c(0,maxAxis))}
    else { p1 <- p1+geom_point(col=colsToUse)+xlim(c(0,maxAxis))}
    p1 <- p1+ylim(c(0,maxAxis))+geom_abline(intercept=0,slope=1, col="red")
    p1 <- p1+xlab("-log10(Expected)")+ylab("-log10(Observed)")
    return(p1)
}

## Create a new sample column based on adding an additional factor
newDat <- function(inDat,
                   inlevel="ctx",
                   colName="Tissue",
                   sampleCol="AnimalID",
                   outCol="SampleID")
{
    inDat[,colName] <- inlevel
    inDat[,outCol] <- paste(as.character(inDat[,sampleCol]),inDat[, colName], sep="_")
    return(inDat)
}

## Pairs plots
panel.cor.scale <- function(x,
                            y,
                            digits=2,
                            prefix="",
                            cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = (cor(x, y,use="pairwise"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * abs(r))
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = (cor(x, y,use="pairwise"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex )
}


panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE)
{
    par(pch ="." )
    if (smooth )
        {
            if (scale)
                {
                    pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth)
                }
            else{
                pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
            } #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
        }
    else #smooth is not true
        {
            if (scale)
                {
                    pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale)
                }
            else
                {
                    pairs(x,diag.panel=panel.hist,upper.panel=panel.cor)
                }
        } #end of else (smooth)
} #end of function

## Get Gene list from threshold
genesToKeep <- function(inDat, threshold)
{
                                        #inDat <- cortexDatFilt1
                                        #threshold <- 12
    inDat <- dcast(gene~SampleID,data = inDat, value.var="value")
    rowsToKeep <- apply(inDat[,-c(1)],MARGIN = 1, FUN =function(x, t) { y <- TRUE ;  if(any(x < t)){ y <- FALSE}; return (y)}, t= threshold)
    genesKept<- as.character(inDat$gene[rowsToKeep])
    return(genesKept)
                                        #rm(inDat)
                                        #rm(threshold)
                                        #rm(rowsToKeep)
                                        #rm(genesKept)
}

### Get alignmentStats data
## column order
## 1  SummStats:Sample
## 3	FR
## 22	RF
## 10	FF
## 18	RR
## 21	numAlignments
## 13	numUniqAlign
## 35	numUniqAlignPE
## 7	numUniqAlignPEProper
## 25	numUniqAlignPEProperFR
## 6	numUniqAlignPEProperFF
## 33	numUniqAlignPEProperRR
## 9	numUniqAlignPELong
## 30	numUniqAlignPELongConc
## 8	numUniqAlignPELongConcFR
## 17	numUniqAlignPELongConcFF
## 36	numUniqAlignPELongConcRR
## 23	numUniqAlignPELongDisc
## 34	numUniqAlignPELongDiscFR
## 28	numUniqAlignPELongDiscFF
## 20	numUniqAlignPELongDiscRR
## 44	numMultipleAlign
## 39	numMultipleAlignPE
## 2	numMultipleAlignPEProper
## 38	numMultipleAlignPEProperFR
## 42	numMultipleAlignPEProperFF
## 5	numMultipleAlignPEProperRR
## 24	numMultipleAlignPELong
## 37	numMultipleAlignPELongConc
## 40	numMultipleAlignPELongConcFR
## 45	numMultipleAlignPELongConcFF
## 19	numMultipleAlignPELongConcRR
## 43	numMultipleAlignPELongDisc
## 4	numMultipleAlignPELongDiscFR
## 11	numMultipleAlignPELongDiscFF
## 26	numMultipleAlignPELongDiscRR
## 16	numUniqAlignSE
## 29	numUnmapped
## 15	numMultipleAlignSE
## 27	numERCCUniqAlignPE
## 12	numERCCUniqAlignDisc
## 14	numERCCUniqAlignPEDisc
## 41	numERCCMultipleAlignPE
## 31	numERCCMultipleAlignDisc
## 32	numERCCMultipleAlignPEDisc

getAlignStats <-function(infile)
{

  dat<- read.csv(infile, header=T)
  dat <- dat [, c(1,3,22,10,18,21,13,35,7,25,6,33,9,30,8,17,36,23,34,28,20,44,39,2,38,42,5,24,37,40,45,19,43,4,11,26,16,29,15,27,12,14,41,31,32)]
  names(dat)[1] <- "SampleID"
  return(dat)
}


## Draw CPM Plots by Genotype
plotbyGtype <- function( inDat, plotFile="plots/7qRegionMeans.png")
{
                                        #inDat <- ctxDat7q.Batch2
                                        #plotFile <- "plots/7qRegionMeansBatch2_%d.png"
    length(unique(inDat$gene))
    print(unique(inDat$gene))
    head(inDat)
    gtypeMeans <- aggregate(cpm~gene:Gtype, data=inDat, mean)
    geneInfo <- inDat[,c("gene","chr","start","end")]
    geneInfo <- geneInfo[!duplicated(geneInfo),]
    gtypeMeans <- merge(gtypeMeans,geneInfo , by="gene")
    gtypeMeans <- gtypeMeans[order(gtypeMeans$start,gtypeMeans$Gtype),]
    head(gtypeMeans)

    gtypeNorm <- subset(gtypeMeans, Gtype=="WT")
    gtypeNorm <- gtypeNorm [ ,c("gene","cpm")]
    names(gtypeNorm)[2] <- "wtcpm"
    print(head(gtypeNorm))

    gtypeMeans <- merge(gtypeMeans, gtypeNorm, by="gene")
    gtypeMeans$log2FC <- log2(gtypeMeans$cpm/gtypeMeans$wtcpm)
    gtypeMeans$FC <- gtypeMeans$cpm/gtypeMeans$wtcpm
    geneNames <- read.csv("genesToNames.csv",header=T)
    head(geneNames)
    gtypeMeans <- merge(gtypeMeans,geneNames, by="gene")

    print(plotFile)

    png(file= plotFile ,height =8, width= 12, res = 600,units ="in" )
    g1 <- ggplot(gtypeMeans,aes(x=as.factor(start/1000),y=log(cpm), colour=Gtype))
    print(g1+geom_point()+theme(axis.text.x=element_text(angle=90)))

    g1 <- ggplot(gtypeMeans,aes(x=as.factor(start/1000),y=log2FC, colour=Gtype))
    print(g1+geom_point()+theme(axis.text.x=element_text(angle=90)))

    g1 <- ggplot(gtypeMeans,aes(x=as.factor(start/1000),y=FC, colour=Gtype))
    print(g1+geom_point()+theme(axis.text.x=element_text(angle=90)))
    dev.off()
}


## Read Alignment Stat datq
readAlignStats <-function(alignFile,nSamp,locDup,locGpRna)
{
    ##alignFile <-"fdmouseSummaries.txt"
    ##locDup <- 125
    ##locGpRna <- 218
    ##nSamp <- 30
    dupStat<- read.csv(alignFile,skip=locDup-1, header=T,nrows=nSamp)
    names(dupStat) <- c("Sample","UNMAPPED_READS" ,"ESTIMATED_LIBRARY_SIZE","OPTICAL_DUPS","readPairs",
                        "UNPAIRED_READS", "DUP_percent","LIBRARY","READ_PAIR_DUPS",
                        "UNPAIRED_DUPS")
    gpRnaSeqQc <-read.csv(alignFile,skip=locGpRna-1, header=T,nrows=nSamp)

    names(gpRnaSeqQc) <-c("SampleID","Total_Reads_Mapped","rRNA","Mapped_Reads", "Unique_Mapping_Rate",
                          "Intragenic_Rate","Note","Mapped_Unique","Intergenic_Rate",
                          "End1_Sense", "Genes_Detected","Num_Gaps","Mean_CV",
                          "End1_Mapping_Rate","End2_Mapping_Rate","Mapped_Pairs","End1_Mismatch_Rate",
                          "Estimated_Library_Size","X5__Norm","Failed_Vendor_QC_Check","End2_Sense",
                          "Mapping_Rate", "Fragment_Length_Mean","Duplication_Rate_of_Mapped",
                          "rRNA_rate","Intronic_Rate","Read_Length", "Transcripts_Detected",
                          "End1_Antisense","Gap_Pct","Expression_Profiling_Efficiency","End1_Pct_Sense",
                          "Unpaired_Reads","Cumul_Pct_Gap_Length","No_Pct_Covered_5Prime",
                          "Alternative_Alignments","Chimeric_Pairs","End2_Mismatch_Rate",
                          "Fragment_Length_StdDev",
                          "Mapped_Unique_Rate_of_Total","Exonic_Rate",
                          "End2_Pct_Sense","Sample",  "Base_Mismatch_Rate","End2_Antisense",
                          "Mean_Per_Base_Cov")

    return(list(dupMetrics=dupStat, gpRnaSeqQc=gpRnaSeqQc))
}

## Read Alignment Stat datq
readGpRNAStats <-function(alignFile)
  {
      gpRnaSeqQc <-read.csv(alignFile,header=T)

      names(gpRnaSeqQc) <-c("SampleID","Total_Reads_Mapped","rRNA","Mapped_Reads", "Unique_Mapping_Rate",
                            "Intragenic_Rate","Note","Mapped_Unique","Intergenic_Rate",
                            "End1_Sense", "Genes_Detected","Num_Gaps","Mean_CV",
                            "End1_Mapping_Rate","End2_Mapping_Rate","Mapped_Pairs","End1_Mismatch_Rate",
                            "Estimated_Library_Size","X5__Norm","Failed_Vendor_QC_Check","End2_Sense",
                            "Mapping_Rate", "Fragment_Length_Mean","Duplication_Rate_of_Mapped",
                            "rRNA_rate","Intronic_Rate","Read_Length", "Transcripts_Detected",
                            "End1_Antisense","Gap_Pct","Expression_Profiling_Efficiency","End1_Pct_Sense",
                            "Unpaired_Reads","Cumul_Pct_Gap_Length","No_Pct_Covered_5Prime",
                            "Alternative_Alignments","Chimeric_Pairs","End2_Mismatch_Rate",
                            "Fragment_Length_StdDev",
                            "Mapped_Unique_Rate_of_Total","Exonic_Rate",
                            "End2_Pct_Sense","Sample",  "Base_Mismatch_Rate","End2_Antisense",
                            "Mean_Per_Base_Cov")

      return(gpRnaSeqQc)
  }
## Gp Rna Setup
setupGpRnaSeq <- function(inDat)
{
    inDat$Mapping_Rate_of_Total <- inDat$Mapped_Reads/(inDat$FastqReads*2)
    inDat$Mapped_Unique_Rate_of_Total <- inDat$Mapped_Unique/(2*inDat$FastqReads)
    inDat$Chimeric_Rate <- inDat$Chimeric_Pairs/inDat$Mapped_Pairs

    mappingCols <- c("SampleID", "FastqReads","Total_Reads_Mapped","Mapped_Reads", "Mapped_Unique", "Mapped_Pairs","Chimeric_Pairs")

    otherMappingCols <-c("SampleID", "Alternative_Alignments", "Unpaired_Reads", "End1_Sense", "End2_Sense", "End1_Antisense",  "End2_Antisense")

    rateCols <- c("SampleID", "Mapping_Rate_of_Total", "Mapping_Rate", "Mapped_Unique_Rate_of_Total",  "Unique_Mapping_Rate","Chimeric_Rate")
    otherRateCols <-c( "SampleID", "End1_Mapping_Rate", "End2_Mapping_Rate", "Base_Mismatch_Rate",
                      "End1_Mismatch_Rate", "End2_Mismatch_Rate", "End1_Pct_Sense", "End2_Pct_Sense")


    otherStatsCols <- c("SampleID","Duplication_Rate_of_Mapped","Mean_Per_Base_Cov","Estimated_Library_Size",
                        "Fragment_Length_Mean","Fragment_Length_StdDev","Genes_Detected", "Transcripts_Detected")

    geneStatsCols<-c("SampleID", "Intragenic_Rate","Exonic_Rate",	"Intronic_Rate","Intergenic_Rate",	"Expression_Profiling_Efficiency")
    return(list(data=inDat,mappingCols=mappingCols,otherMappingCols=otherMappingCols,rateCols=rateCols,otherRateCols=otherRateCols,
                otherStatCols=otherStatsCols,geneStatsCols=geneStatsCols))

}

## MySQL functions

uploadCountsToMySQLDB <- function(inCountsFile, strand="sense", type="counts", mapping="multiple",tool="bedtools",desc="NULL" )
{
    library(reshape2)
    library(RMySQL)
                                        #inCountsFile <- "/Users/ashok/Documents/Research/CHGR/16pMusTissueSpecific/DataAnalysis/Counts/bfatGeneCounts_sense.sum.cnts.csv"
                                        #type <-"counts"
                                        #strand <- "sense"
                                        #mapping <- "multiple"
                                        #tool <- "bedtools"
                                        #desc <- "NULL"
    counts_dat <- read.csv(inCountsFile, header=T)
    counts_dat_long <- melt(counts_dat,value.name = "counts",id.vars = c("gene","chr","start","end"))
    counts_dat_long$location=paste(counts_dat_long$chr, paste(counts_dat_long$start,counts_dat_long$end, sep="-"), sep=":")
    counts_dat_long <- subset(counts_dat_long, select=c("gene","variable","counts","location"))
    names(counts_dat_long) <-c("geneid","id","expression","location")
    counts_dat_long$type<-type
    counts_dat_long$strand <- strand
    counts_dat_long$mapping <- mapping
    counts_dat_long$tool <- tool
    counts_dat_long$desc <- desc
    print(head(counts_dat_long))
    con <- dbConnect(drv=MySQL(), user='talkLabAdmin', password ='_Talk0w_',dbname='TalkowskiGenomicsDB', host='mysql2.dipr.partners.org')
    dbWriteTable(conn = con, name = "geneRNASeq",value = counts_dat_long,append=T, row.names=F)
    dbDisconnect(con)
}


## PCA plot function
## This function is used to plot a PC object selecting specific columns and a color variable


myPlotPC <- function(pr.obj,col1,col2, inCol=c(rep("red",4),rep("blue",4)))
    {
        plot(pr.obj$x[,col1]/pr.obj$sdev[col1],pr.obj$x[,col2]/pr.obj$sdev[col2],
             pch=19, col=inCol,
             xlab=paste("PC",col1,sep=''),ylab=paste("PC",col2,sep=''))
        text(pr.obj$x[,col1]/pr.obj$sdev[col1],pr.obj$x[,col2]/pr.obj$sdev[col2],
             labels=rownames(pr.obj$x),cex=0.7,
             col=inCol, pos=1)
    }


## ********************************
## This Section is RNASeq Specific

## Create Manhattan plot setup

## This function takes as input an ERCC genes list data in Long format and merges it with the ERCC
## Concentrations etc for each ercc gene. Two important considerations are:
##  1. input Data should have the gene identifier column with the geneID
##  2. It also creates a merged Ercc concentration colum for downstream analysis

mergedErccInfo<- function(inDat)
    {
        erccConcDat <- read.table("~/Documents/Research/CHGR/RNoteBooks/ERCCSpikeIn.txt", header=T)
        head(erccConcDat)
        erccDat <- merge(inDat, erccConcDat, by.x="gene", by.y="ERCCID")
        erccDat$concMix <- erccDat$concMix1
        erccDat$concMix[which(erccDat$Mix == "Mix2")] <- erccDat$concMix2[which(erccDat$Mix == "Mix2")]
        head(erccDat)
        return(erccDat)
    }


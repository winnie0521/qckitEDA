## Create a list of dtracks for plotting with Gviz
makeTracksToPlot <- function(inDat, geneModels,inGenome="mm10",chr=4, inStrand = "+", geneName = "Ikbap endo" , plot=FALSE,setylim=F, ylimits=c(0,200))
{
  dtrack <- list()
  for( i in 1 : dim(bamFiles)[1])
  {
    if (!setylim)
    {
    dtrack[i] <- DataTrack(range=inDat$file[i] ,type="histogram", 
                           name=inDat$sample[i], genome = "mm10", chromosome=chr, strand=inStrand,
                             options(ucscChromosomeNames=FALSE))    
    }else
    {
      dtrack[i] <- DataTrack(range=inDat$file[i] ,type="histogram", 
                                  name=inDat$sample[i], genome = "mm10", chromosome=chr, strand=inStrand,
                                  options(ucscChromosomeNames=FALSE),ylim=ylimits)          
    }
  }
  dtrack[i+1] <- GenomeAxisTrack()
  dtrack[i+2] <- GeneRegionTrack(geneModels, genome=inGenome,chromosome=chr,name= geneName,
                                 options(ucscChromosomeNames=FALSE))
  if(plot)
  {plotTracks(dtrack, showId=T)}
  return(dtrack)
}
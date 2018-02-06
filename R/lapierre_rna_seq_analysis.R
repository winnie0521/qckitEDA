library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gplots)
library(RSQLite)
library(dendextend)
library(ascii)
source('R/generalfunctions.R', echo=TRUE)

dat <- read.csv("Data/counts.csv", header=T)
rows_to_keep= apply(dat[,-c(1)],1,FUN=function(x)
      {
          y <- 0;
          if (all(x > 0 )) y <-1;
          return(y)
      }
      )

dat_filt <- list()
dat_filt[["wide"]]<- dat[as.logical(rows_to_keep > 0),]
dat_filt[["long"]]<- melt(dat_filt$wide, id="gene", variable.name="Sample")

print (dim(dat)[1], dim(dat_filt)[1])

## Load Design
model.mat <- read.csv("Data/expt_design.csv", header=T)
row.names(model.mat) <- model.mat$Sample
dat_filt[["long"]]<- merge(dat_filt$long, model.matrix, by="Sample", all=T)

## Convert Counts to CPM
libSizes <- aggregate(data=dat_filt$long, value~Sample,FUN=sum)
names(libSizes)[2] <- "libsize"
dat_filt$long <- merge(dat_filt$long, libSizes, by="Sample", all=T)
dat_filt$long$cpm <- dat_filt$long$value*1000000/dat_filt$long$libsize

#Add Gene annotations
library(biomaRt)
ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")
filters = listFilters(ensembl)
attrs.to.get = listAttributes(ensembl)
common_ids = getBM(attributes="ensembl_gene_id", filters='wormbase_gene',value=dat_filt$wide$gene, mart = ensembl)
annoDat <- getBM(attributes=c("ensembl_gene_id",attrs.to.get[c(15,6:9,20,21,5),1]),filters = 'ensembl_gene_id',values =common_ids, mart = ensembl)
row.names(annoDat) <- annoDat$ensembl_gene_id

g_1 <- ggplot(annoDat, aes(x=gene_biotype))+geom_bar(stat="count")
g_1<- ggplotly(g_1)
htmlwidgets::saveWidget(g_1,"plots/lapierre_plotly_gene_biotype.html")

g1 <- ggplot(dat_filt$long,aes(x=Sample,y=cpm, fill=Gtype))+geom_boxplot()+scale_y_log10()
g1 <- g1+theme(axis.text.x=element_text(size=8, angle=90))+xlab("")
print(g1)

print(g1)

g1 <- ggplot(dat_filt$long,aes(x=Sample,y=log(cpm+0.0000001), fill=Gtype))+geom_boxplot()
g1+theme(axis.text.x=element_text(size=8, angle=90))+xlab("")
print(g1)

print(g1)

tmpDat.cpm <- dcast(data=dat_filt$long,gene~Sample, value.var ="cpm")
g1<- ggMDSplot(tmpDat.cpm,modMat =model.matrix, modCol = 2,sampleLoc = "Sample", txtSize = 3)
print(g1)

print(g1)
rm(g1)
rm(tmpDat.cpm)

dat_filt$long$logcpm <- log(dat_filt$long$cpm+0.00001)
tmpDat.cpm <- dcast(data=dat_filt$long,gene~Sample, value.var ="logcpm")
g1<- ggMDSplot(tmpDat.cpm,modMat =model.matrix, modCol = 2,sampleLoc = "Sample", txtSize = 3)
print(g1)

print(g1)
rm(g1)
rm(tmpDat.cpm)

tmpDat.cpm <- dcast(dat_filt$long,gene~Sample,value.var = "cpm")[,-c(1)]
geneCor.pr <- cor(tmpDat.cpm, use="complete.obs", method="pearson")
geneCor.sp <- cor(tmpDat.cpm, use="complete.obs", method="spearman")
geneCor.log.pr <- cor(log(tmpDat.cpm+0.000001), use="complete.obs", method="pearson")
geneCor.log.sp <- cor(log(tmpDat.cpm+0.000001), use="complete.obs", method="spearman")
my_palette <- colorRampPalette(rev(brewer.pal(9,"RdGy")))(n = 300)

heatmap.2(geneCor.pr, trace="none", dendrogram="column",  col=my_palette, scale="none")

heatmap.2(geneCor.pr, trace="none", dendrogram="column",  col=my_palette, scale="none")

heatmap.2(geneCor.sp, trace="none", dendrogram="column", col=my_palette, scale="none")

heatmap.2(geneCor.sp, trace="none", dendrogram="column", col=my_palette, scale="none")

heatmap.2(geneCor.log.pr, trace="none", dendrogram="column",  col=my_palette, scale="none")

heatmap.2(geneCor.log.pr, trace="none", dendrogram="column",  col=my_palette, scale="none")

heatmap.2(geneCor.log.sp, trace="none", dendrogram="column", col=my_palette, scale="none")

heatmap.2(geneCor.log.sp, trace="none", dendrogram="column", col=my_palette, scale="none")
rm(tmpDat.cpm)

tmpDat.cpm.log <- dcast(dat_filt$long,gene~Sample,value.var = "logcpm")[,-c(1)]
names(tmpDat.cpm.log)
pairs.panels(tmpDat.cpm.log)

pairs.panels(tmpDat.cpm.log)

tmpDat.cpm <- dcast(dat_filt$long,gene~Sample,value.var = "cpm")[,-c(1)]
tmpDat.cpm <- t(tmpDat.cpm)
tmpDat.cpm <- as.matrix(tmpDat.cpm)
colnames(tmpDat.cpm)
tmpDat.cpm <- dist(tmpDat.cpm)
par(bg="lightyellow")
dend <- as.dendrogram(hclust(tmpDat.cpm,method="average"))
##plot(hclust(tmpDat.cpm,method="average"))
colors_to_use <- subset(dat_filt$long, gene==levels(factor(dat_filt$long$gene))[1])$Gtype
colors_to_use <- gsub("WT","black", colors_to_use)
colors_to_use <- gsub("glp_1_vit_2_GFP","blue", colors_to_use)
colors_to_use <- gsub("vit_2_GFP","green", colors_to_use)
colors_to_use <- gsub("glp_1","red", colors_to_use)
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use
plot(dend)

plot(dend)

par(bg="lightyellow")
##plot(hclust(tmpDat.cpm,method="ward.D"))
dend <- as.dendrogram(hclust(tmpDat.cpm,method="ward"))
##plot(hclust(tmpDat.cpm,method="average"))
colors_to_use <- subset(dat_filt$long, gene==levels(factor(dat_filt$long$gene))[1])$Gtype
colors_to_use <- gsub("WT","black", colors_to_use)
colors_to_use <- gsub("glp_1_vit_2_GFP","blue", colors_to_use)
colors_to_use <- gsub("vit_2_GFP","green", colors_to_use)
colors_to_use <- gsub("glp_1","red", colors_to_use)
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use
plot(dend)

plot(dend)
rm(tmpDat.cpm)

tmpDat.cpm <- dcast(dat_filt$long,gene~Sample,value.var = "logcpm")[,-c(1)]
tmpDat.cpm <- t(tmpDat.cpm)
tmpDat.cpm <- as.matrix(tmpDat.cpm)
colnames(tmpDat.cpm)
tmpDat.cpm <- dist(tmpDat.cpm)
par(bg="lightyellow")
dend <- as.dendrogram(hclust(tmpDat.cpm,method="average"))
##plot(hclust(tmpDat.cpm,method="average"))
colors_to_use <- subset(dat_filt$long, gene==levels(factor(dat_filt$long$gene))[1])$Gtype
colors_to_use <- gsub("WT","black", colors_to_use)
colors_to_use <- gsub("glp_1_vit_2_GFP","blue", colors_to_use)
colors_to_use <- gsub("vit_2_GFP","green", colors_to_use)
colors_to_use <- gsub("glp_1","red", colors_to_use)
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use
plot(dend)

plot(dend)

par(bg="lightyellow")
##plot(hclust(tmpDat.cpm,method="ward.D"))
dend <- as.dendrogram(hclust(tmpDat.cpm,method="ward"))
##plot(hclust(tmpDat.cpm,method="average"))
colors_to_use <- subset(dat_filt$long, gene==levels(factor(dat_filt$long$gene))[1])$Gtype
colors_to_use <- gsub("WT","black", colors_to_use)
colors_to_use <- gsub("glp_1_vit_2_GFP","blue", colors_to_use)
colors_to_use <- gsub("vit_2_GFP","green", colors_to_use)
colors_to_use <- gsub("glp_1","red", colors_to_use)
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use
plot(dend)

plot(dend)
rm(tmpDat.cpm)

## Load Libraries
library(DSS)
library(edgeR)
library(biomaRt)

## Setup the data
X<-model.matrix(~Gtype, data=mod.mat)
## For now column order matches the model.matrix
dat_filt[["matrix"]] <- as.matrix(dat_filt$wide[, -c(1)])
rownames(dat_filt$matrix) <- dat_filt$wide$gene
head(dat_filt$matrix)
seqData <- newSeqCountSet(dat_filt$matrix, as.data.frame(X))
seqData <- estNormFactors(seqData)
seqData <- estDispersion(seqData)

##fit GLM
fit.edgeR <- glmFit(dat_filt$matrix, X, lib.size=normalizationFactor(seqData),dispersion=dispersion(seqData))
fitql.edgeR <- glmQLFit(dat_filt$matrix, X, lib.size=normalizationFactor(seqData),dispersion=dispersion(seqData))

pair_vector = sprintf("%s-%s", "Intercept", "GtypeWT") # Samples to be compared
pair_contrast = makeContrasts(contrasts=pair_vector, levels=X) # Make contrast

edgeR.F <-  glmQLFTest(glmfit=fitql.edgeR, contrast=pair_contrast)
resTab <- edgeR.F$table
resTab$wormbase_gene <- rownames(edgeR.F$table)
resTab<-merge(resTab,annoDat,by="wormbase_gene", all.x=T)

g1<-qqGGplot(resTab$PValue)
#g1 <-ggplotly(g1, tooltip=resTab$external_gene_name[order(resTab$PValue, decreasing=F)])
print(g1)

g1 <- ggplot(resTab, aes(y=-log10(PValue),x=logFC, text=external_gene_name))
g1<- g1+ geom_point(alpha=1/2,aes(colour = cut(PValue, c(0, 0.001, 0.05, 1))))
g1 <- g1+ scale_color_manual(name = "pvalue",
                     values = c("(0,0.001]" = "red",
                                  "(0.001,0.05]" = "yellow",
                                  "(0.05, 1]" = "blue"),
                     labels = c("<= 0.001", "0.001 < pvalue <= 0.05", "> 0.05"))
print(g1)

g11<-ggplotly(g1)
htmlwidgets::saveWidget(g11,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_plotly_glp1_wt_volcano.html")

## tableHTML approach
## t1<-tableHTML(resTab, theme="scientific")
## write_tableHTML(t1, file="lapierre_rna_seq_analysis_tmptable1.html")

## system2("sed", args=c('s/class=table_resTab/id=\\"glp1_wt\\"/',"lapierre_rna_seq_analysis_tmptable1.html"," > lapierre_rna_seq_analysis_tmptable11.html"))
## system2("sed", args=c( "'s/border=0>/class=\"display\" cellspacing=\"0\" width=\"100%\">/'","lapierre_rna_seq_analysis_tmptable11.html"," > lapierre_rna_seq_analysis_table1.html"))

tab1 <- resTab[order(resTab$PValue),]
tab1 <- datatable(tab1[1:1000,])
htmlwidgets::saveWidget(tab1,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_dt_glp1_wt_table.html")

pair_vector = sprintf("%s-%s", "Gtypevit_2_GFP", "GtypeWT") # Samples to be compared
pair_contrast = makeContrasts(contrasts=pair_vector, levels=X) # Make contrast

edgeR.F <-  glmQLFTest(glmfit=fitql.edgeR, contrast=pair_contrast)
resTab2 <- edgeR.F$table
resTab2$wormbase_gene <- rownames(edgeR.F$table)
resTab2<-merge(resTab2,annoDat,by="wormbase_gene", all.x=T)

g2<-qqGGplot(resTab2$PValue)
print(g2)

g21 <- ggplot(resTab2, aes(y=-log10(PValue),x=logFC, text=external_gene_name))+geom_point()
g21<- g21+ geom_point(alpha=1/2,aes(colour = cut(PValue, c(0, 0.001, 0.05, 1))))
g21 <- g21+ scale_color_manual(name = "pvalue",
                     values = c("(0,0.001]" = "red",
                                  "(0.001,0.05]" = "yellow",
                                  "(0.05, 1]" = "blue"),
                     labels = c("<= 0.001", "0.001 < pvalue <= 0.05", "> 0.05"))
print(g21)

g211<-ggplotly(g21)
htmlwidgets::saveWidget(g211,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_plotly_vit2_wt_volcano.html")

tab2 <- resTab2[order(resTab2$PValue),]
tab2 <- datatable(tab2[1:1000,])
htmlwidgets::saveWidget(tab2,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_dt_vit2_wt_table.html")

pair_vector = sprintf("%s-%s", "Gtypevit_2_GFP", "Gtypeglp_1_vit_2_GFP") # Samples to be compared
pair_contrast = makeContrasts(contrasts=pair_vector, levels=X) # Make contrast

edgeR.F <-  glmQLFTest(glmfit=fitql.edgeR, contrast=pair_contrast)
resTab3 <- edgeR.F$table
resTab3$wormbase_gene <- rownames(edgeR.F$table)
resTab3<-merge(resTab3,annoDat,by="wormbase_gene", all.x=T)

g3<-qqGGplot(resTab3$PValue)
#g1 <-ggplotly(g1, tooltip=resTab$external_gene_name[order(resTab$PValue, decreasing=F)])
print(g3)

g31 <- ggplot(resTab3, aes(y=-log10(PValue),x=logFC, text=external_gene_name))+geom_point()
g31<- g31+ geom_point(alpha=1/2,aes(colour = cut(PValue, c(0, 0.001, 0.05, 1))))
g31 <- g31+ scale_color_manual(name = "pvalue",
                     values = c("(0,0.001]" = "red",
                                  "(0.001,0.05]" = "yellow",
                                  "(0.05, 1]" = "blue"),
                     labels = c("<= 0.001", "0.001 < pvalue <= 0.05", "> 0.05"))
print(g31)

g311<-ggplotly(g31)
htmlwidgets::saveWidget(g311,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_plotly_vit2_vit2glp1_volcano.html")

tab3 <- resTab3[order(resTab3$PValue),]
tab3 <- datatable(tab3[1:1000,])
htmlwidgets::saveWidget(tab3,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_dt_vit2_vit2glp1_table.html")

pair_vector = sprintf("%s-%s", "Intercept", "Gtypeglp_1_vit_2_GFP") # Samples to be compared
pair_contrast = makeContrasts(contrasts=pair_vector, levels=X) # Make contrast

edgeR.F <-  glmQLFTest(glmfit=fitql.edgeR, contrast=pair_contrast)
resTab4 <- edgeR.F$table
resTab4$wormbase_gene <- rownames(edgeR.F$table)
resTab4<-merge(resTab4,annoDat,by="wormbase_gene", all.x=T)

g4<-qqGGplot(resTab4$PValue)
#g1 <-ggplotly(g1, tooltip=resTab$external_gene_name[order(resTab$PValue, decreasing=F)])
print(g4)

g41 <- ggplot(resTab4, aes(y=-log10(PValue),x=logFC, text=external_gene_name))+geom_point()
g41<- g41+ geom_point(alpha=1/2,aes(colour = cut(PValue, c(0, 0.001, 0.05, 1))))
g41 <- g41+ scale_color_manual(name = "pvalue",
                     values = c("(0,0.001]" = "red",
                                  "(0.001,0.05]" = "yellow",
                                  "(0.05, 1]" = "blue"),
                     labels = c("<= 0.001", "0.001 < pvalue <= 0.05", "> 0.05"))
print(g41)

g411<-ggplotly(g41)
htmlwidgets::saveWidget(g411,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_plotly_glp1_vit2glp1_volcano.html")

tab4 <- resTab4[order(resTab4$PValue),]
tab4 <- datatable(tab4[1:1000,])
htmlwidgets::saveWidget(tab4,"/Users/aragaven/Documents/Research/RNotes/plots/lapierre_dt_glp1_vit2glp1_table.html")

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("edgeR","limma","calibrate"))
#biocLite(c("Glimma"))
library(edgeR)
library(gplots) #install.packages("gplots")
library(RColorBrewer)
library(limma)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/POP#2/RNA/DEG/")

######################## Filtered counts file
Counts <- read.csv("1.Raw_counts/20180321_Counts_POP2.csv",
                   header = T,stringsAsFactors = F,row.names = 1)
dim(Counts)

######################## Making the metadata table
## This will only work if your sample names are separated by an underscore (_)
meta <- data.frame(do.call("rbind",strsplit(colnames(Counts),"_")),
                        row.names = colnames(Counts),stringsAsFactors = F)
## Manually set the column names
colnames(meta) <- c("Genotype","Time","Soil","Treatment","R")
head(meta)
### Important: the rownames of the metadata table should match the column names of your counts table


####################### Sanity check of counts and meta data compatibility
checkNames <- function(){
  if ((all(colnames(Counts) %in% rownames(meta)))){
    cat ("all is good, # of samples:", nrow(meta))
  } else {
    warning("Check sample names")
  }}
###
checkNames()

####################### Filter individual datasets

## Filter by Genotype
genoFilter <- T
if (genoFilter){
  cat ("Filtering by genotype \n")
  idxGeno <- grep("SQR",meta$Genotype)
  Counts <- Counts[,idxGeno]
  meta <- meta[idxGeno,]
} else { cat ("No filter \n") }
checkNames()


## Filter by Day
timeFilter <- T
if (timeFilter){
  cat ("Filtering by genotype \n")
  idxTime <- grep("2",meta$Time)
  Counts <- Counts[,idxTime]
  meta <- meta[idxTime,]
} else { cat ("No filter \n") }
checkNames()

####################### Removing lowly-expressed genes
dge <- DGEList(counts=Counts,remove.zeros = T)
sampleMin <- 3 
minCPM <- 1
isexpr <- rowSums(cpm(dge) > minCPM) >= sampleMin 
dge <- dge[isexpr,,keep.lib.size = FALSE]

#######################  Calculate the normalization factors with TMM
dge.norm <- calcNormFactors(dge, method = "TMM")
dge.norm$samples$norm.factors

par(mfrow=c(1,2))
lcpm <- cpm(dge, log=TRUE)
boxplot(lcpm)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm2 <- cpm(dge.norm, log=TRUE)
boxplot(lcpm2)
title(main="B. Example: Normalised data",ylab="Log-cpm")

####################### Differential expression - design matrix
meta$Soil <- as.factor(meta$Soil)
meta$Soil <- relevel(meta$Soil,ref = "S")
meta$Treatment <- as.factor(meta$Treatment)
meta$Treatment <- relevel(meta$Treatment,ref = "CTR")

design <- model.matrix(~Soil*Treatment, data=meta)
colnames(design) <- gsub("Soil|Day|Name|Genotype|Treatment|:|-|/|Groups","",colnames(design))
head(design)

##################### MDS plots
meta$ST <- paste(meta$Soil,meta$Treatment,sep="_")
meta$ST <- as.factor(meta$ST)

par(mfrow=c(1,3))
col.soil <- meta$Soil
levels(col.soil) <-  brewer.pal(nlevels(col.soil), "Set1")
col.soil <- as.character(col.soil)
col.tr <- meta$Treatment
levels(col.tr) <-  brewer.pal(nlevels(col.tr), "Set2")
col.tr <- as.character(col.tr)
col.sampl <- meta$ST
levels(col.sampl) <-  brewer.pal(nlevels(col.sampl), "Set3")
col.sampl <- as.character(col.sampl)

plotMDS(lcpm2, labels=meta$Soil,col=col.soil)
title(main="A. Soil")
plotMDS(lcpm2, labels=meta$Treatment, col=col.tr, dim=c(3,4))
title(main="B. Treatment")
plotMDS(lcpm2, labels=meta$ST,col=col.sampl)
title(main="C. Sample")

####################### Differential expression - transform with voom (and renormalize with quantile).
v <- voomWithQualityWeights(dge.norm, design, plot = TRUE, normalize.method = "quantile")
# v <- voom(dge.norm, design, plot = TRUE,normalize.method = "quantile")

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

####################### Save Expression Data
cpmExpression <- cpm(dge.norm)
voomExpression <- v$E
colnames(cpmExpression) <- rownames(meta)
colnames(voomExpression) <-  rownames(meta)

save(cpmExpression,file = "0.DEG2021/SQR2/cpmExpression_SQR2.RData")
save(voomExpression,file = "0.DEG2021/SQR2/voomExpression_SQR2.RData")
write.table(cpmExpression, file = "0.DEG2021/SQR2/voomExpression_SQR2.csv",sep = ",",quote = F,col.names = T )


##########################  DGE
v2 <- lmFit(v, design)
fit <- eBayes(v2)
colSums(abs(decideTests(fit)))

coefficientsTested <- colnames(fit$coefficients)[-1]

DEList <- lapply(coefficientsTested,function(x){
  DEresults <- topTable(fit, coef=x, adjust="BH",number = Inf,sort.by = "none")
  colnames(DEresults) <- paste(x,colnames(DEresults),sep = ".")
  return(DEresults)
})
names(DEList) <- coefficientsTested
names(DEList)

## Get only significant genes

sapply(DEList,function(x){
  nrow(x[x[,grep("adj.P.Val",colnames(x))] < 0.05,])
  #nrow(x[x[,grep("qValue",colnames(x))] < 0.05,])
})


significantDE <- lapply(DEList,function(x){ x[x[,grep("adj.P.Val",colnames(x))] < 0.05,] })
sapply(significantDE,nrow)
shortName <- "SQR2"
saveDElist <- paste0("DEList_",shortName,".RData")
save(file = saveDElist,DEList)
sapply(DEList,colnames)

names(significantDE)

write.table( significantDE[["N"]],file = "0.DEG2021/SQR2/SQR2_Soil.csv",sep = ",",quote = F,col.names = T)
write.table( significantDE[["STRIGA"]],file = "0.DEG2021/SQR2/SQR2_Striga.csv",sep = ",",quote = F,col.names = T)
write.table( significantDE[["NSTRIGA"]],file = "0.DEG2021/SQR2/SQR2_SoilXStriga.csv",sep = ",",quote = F,col.names = T)

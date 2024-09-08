##################################################################################
## Figure 2(A): 
## Create network modules, visualize the modules, associate with anthropometry
## !!!!!MM and GS analysis were not weighted!!!!!!!
##################################################################################
#install.packages("WGCNA")
#BiocManager::install(c( "AnnotationDbi", "impute", "GO.db","preprocessCore"))
#NOTE: fixing heatmap color inversion: https://www.biostars.org/p/394615/
library(WGCNA)
library(impute)
library(nlme)
library(cluster)
library(qvalue)
library(flashClust)
library(Hmisc)
library(dynamicTreeCut)
setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")

rm(list = ls())
anno <- read.csv("../Data/Lipidomics_annotation_file.csv",header=T, sep=",")
anno$met_ID <- chartr(old = ".", new = "_", anno$lipid_identifier)
met_index <- c(4:331)
lipid <- read.csv("../Data/lipid_imputed_wosmoke.csv")
#for (visit in c(0,1,2,4)) {
  
  DF <- lipid[lipid$Visit==0,]
  DFmet <- DF[ ,met_index]
  
  #Choose a set of soft-thresholding powers
  softPower = 10
  # for(k in 1:3) assign(paste0("rank",k), rank(colMeans(DF[DF$batch==k,met_index])))
  # verboseScatterplot(rank1,rank2)
  # verboseScatterplot(rank2,rank3)
  # verboseScatterplot(rank1,rank3)
  for(k in 1:3) assign(paste0("adjacency",k), adjacency(DF[DF$batch==k,met_index], power = softPower, type = "signed"))
  diag(adjacency1)  <- diag(adjacency2) <- diag(adjacency3) <- 0
  for(k in 1:3) assign(paste0("dissTOM",k), 1-TOMsimilarity(get(paste0("adjacency",k)), TOMType="signed"))
  for(k in 1:3) assign(paste0("lipTree",k), flashClust(as.dist(get(paste0("dissTOM",k))), method="average"))
  
  
  par(mfrow=c(1,3))
  plot(lipTree1,xlab="",sub="",main="Metabolites clustering (Batch 1)",
       labels=FALSE,hang=0.04)
  plot(lipTree2,xlab="",sub="",main="Metabolites clustering (Batch 2)",
       labels=FALSE,hang=0.04)
  plot(lipTree3,xlab="",sub="",main="Metabolites clustering (Batch 3)",
       labels=FALSE,hang=0.04)
  
  mColorh = NULL
  for(ds in 0:3){
    tree = cutreeHybrid(dendro = lipTree2, pamStage=FALSE,
                        minClusterSize = (15-3*ds), cutHeight = 0.99,
                        deepSplit = ds, distM = dissTOM2)
    mColorh=cbind(mColorh,labels2colors(tree$labels));
  }
  
  plotDendroAndColors(lipTree2, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE); 
  
  module2 <- mColorh[,1]
  PCs2 = moduleEigengenes(DF[DF$batch==2,met_index], colors=module2)
  ME_2 = PCs2$eigengenes
  distPC2 = 1-abs(cor(ME_2,use="p"))
  distPC2 = ifelse(is.na(distPC2), 0, distPC2)
  pcTree2 = hclust(as.dist(distPC2),method="a")
  MDS_2 = cmdscale(as.dist(distPC2),2)
  colors2 = names(table(module2)) 
  
  par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
  plot(pcTree2, xlab="",ylab="",main="",sub="")
  plot(MDS_2, col= colors2, main="MDS plot", cex=2, pch=19)
  # orderlip = lipTree2$order
  # plotMat(scale(log(DFmet[DF$batch==2,orderlip])) , rlabels= module2[orderlip], clabels=
  #            colnames(DFmet), rcols=module2[orderlip]) 
  # for (which.module in names(table(module2))){
  #   ME = ME_2[, paste("ME",which.module, sep="")]
  #   barplot(ME, col=which.module, main="", cex.main=2,
  #           ylab="eigengene expression",xlab="array sample")
  # }; 
  
  # find the preservation part for other two batches
  par(mfrow=c(1,2))
  plotDendroAndColors(lipTree1, module2, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                      guideHang=0.05, main="Dendrogram and module colors (Batch1)")
  plotDendroAndColors(lipTree2, module2, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                      guideHang=0.05, main="Dendrogram and module colors (Batch2)") 
  
  multiExpr = list(A1=list(data=DFmet[DF$batch==2,]),A2=list(data=DFmet[DF$batch==3,]))
  multiColor = list(A1 = module2)
  mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                        nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400)
  stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
  stats[order(-stats[,2]),c(1:2)] 
  
  
  Cov <- data.frame(ID = DF$ID, pairID = DF$PairID, v = (1/DF$wt)^2, GA_blood = DF$GA_blood,
                    age = DF$age, race = DF$race, preBMI = DF$preBMI, parity = DF$parity, 
                    female = DF$female, GADel = DF$GADel, glucose = DF$Glucose, anthro = DF$bw/1000)
  
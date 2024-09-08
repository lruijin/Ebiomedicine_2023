library(WGCNA)
library(impute)
library(nlme)
library(cluster)
library(ggplot2)
library(dplyr)
rm(list = ls())
setwd("C:/Users/lur5.NIH/Documents/Project4_Lipidomic")
lipid <- read.csv("Data/lipid_imputed_wosmoke_combat.csv", header = T)
anno <- read.csv("Data/Lipidomics_annotation_complete.csv",header=T, sep=",")

met_index <- 4:331

####--------------- Step 1: data cleaning and outliler removal ---------------####
# Prepare settings:
options(stringsAsFactors = FALSE)
lipid0 = lipid[lipid$Visit == 0,]
lipid1 = lipid[lipid$Visit == 1,]
lipid2 = lipid[lipid$Visit == 2,]
lipid4 = lipid[lipid$Visit == 4,]

# We work with four sets:
nSets = 4;

# For easier labeling of plots, create a vector holding descriptive names of the four sets.
setLabels = c("Visit 0", "Visit 1", "Visit 2", "Visit 4")
shortLabels = c("V0", "V1", "V2", "V4")

#Form multi-set expression data: columns corresponding to met_index conain actual expression data.
multiExpr = vector(mode="list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(lipid0[,met_index]))
rownames(multiExpr[[1]]$data) = lipid0$ID;
names(multiExpr[[1]]$data) = names(lipid0)[met_index]

multiExpr[[2]] = list(data = as.data.frame(lipid1[,met_index]))
rownames(multiExpr[[2]]$data) = lipid1$ID;
names(multiExpr[[2]]$data) = names(lipid1)[met_index]

multiExpr[[3]] = list(data = as.data.frame(lipid2[,met_index]))
rownames(multiExpr[[3]]$data) = lipid2$ID;
names(multiExpr[[3]]$data) = names(lipid2)[met_index]

multiExpr[[4]] = list(data = as.data.frame(lipid4[,met_index]))
rownames(multiExpr[[4]]$data) = lipid4$ID;
names(multiExpr[[4]]$data) = names(lipid4)[met_index]

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize

# check that all genes and samples have sufficiently low numbers of missing values"
gag = goodSamplesGenesMS(multiExpr, verbose = 3)
gag$allOK

#cluster the samples on their Euclidean distance, separtely in each set.
sampleTrees = list()
for(set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

#view the four dendrograms:
pdf(file = "Figures/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off()

# The only outlier observed is the one in visit 0. Then we set the base cut line at 22.5 to cut it out.
baseHeight = 19
# Adjust the cut height for the other three datasets for the number of samples
cutHeights = baseHeight * exprSize$nSamples[c(1:4)]/exprSize$nSamples[4]
# Re-plot the dendrograms including the cut lines
pdf(file = "Figures/SampleClustering.pdf", width = 12, height = 12)
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h = cutHeights[set], col = "red")
}
dev.off()

# the actual outlier removal:
Keep_index = vector(mode = "list", length = nSets)
for(set in 1:nSets){
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels == 1)
  #keep = T
  Keep_index[[set]] = keep
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep,]
}
collectGarbage()
exprSize = checkSets(multiExpr)
exprSize

# set the traits data
Traits = vector(mode = "list", length = nSets)
DF = lipid0[Keep_index[[1]],]
Traits[[1]] = list(data = data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                                     age = DF$age, race = factor(as.integer(DF$race == 3)), education = factor(as.integer(DF$Education > 3), ordered = T), 
                                     preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                                     bwz = DF$bwz, bw = DF$bw.x, Length = DF$Length, SSF = DF$SSF, ab_circ = DF$ab_circ,
                                     head_circ = DF$head_circ, diet = DF$diet, 
                                     GWGv0 = DF$GWGv0, GWGv1 = DF$GWGv1, GWGv2 = DF$GWGv2, GWGv4 = DF$GWGv4,
                                     LBW = factor(as.numeric(DF$LBW)), SGA = factor(DF$SGA), LGA = factor(DF$LGA), macro = factor(as.numeric(DF$macro))))
DF = lipid1[Keep_index[[2]],]
Traits[[2]] = list(data = data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                                     age = DF$age, race = factor(as.integer(DF$race == 3)), education = factor(as.integer(DF$Education > 3), ordered = T), 
                                     preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                                     bwz = DF$bwz, bw = DF$bw.x, Length = DF$Length, SSF = DF$SSF, ab_circ = DF$ab_circ,
                                     head_circ = DF$head_circ, diet = DF$diet,
                                     LBW = factor(as.numeric(DF$LBW)), SGA = factor(DF$SGA), LGA = factor(DF$LGA), macro = factor(as.numeric(DF$macro))))
DF = lipid2[Keep_index[[3]],]
Traits[[3]] = list(data = data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                                     age = DF$age, race = factor(as.integer(DF$race == 3)), education = factor(as.integer(DF$Education > 3), ordered = T), 
                                     preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                                     bwz = DF$bwz, bw = DF$bw.x, Length = DF$Length, SSF = DF$SSF, ab_circ = DF$ab_circ,
                                     head_circ = DF$head_circ, diet = DF$diet,
                                     LBW = factor(as.numeric(DF$LBW)), SGA = factor(DF$SGA), LGA = factor(DF$LGA), macro = factor(as.numeric(DF$macro))))
DF = lipid4[Keep_index[[4]],]
Traits[[4]] = list(data = data.frame(ID = DF$ID, pairID = DF$PairID, v = 1/DF$wt, GA_blood = DF$GA_blood,
                                     age = DF$age, race = factor(as.integer(DF$race == 3)), education = factor(as.integer(DF$Education > 3), ordered = T), 
                                     preBMI = DF$preBMI, nulli = factor(as.integer(DF$Nulliparious)), Gender = factor(DF$Gender),
                                     bwz = DF$bwz, bw = DF$bw.x, Length = DF$Length, SSF = DF$SSF, ab_circ = DF$ab_circ,
                                     head_circ = DF$head_circ, diet = DF$diet,
                                     LBW = factor(as.numeric(DF$LBW)), SGA = factor(DF$SGA), LGA = factor(DF$LGA), macro = factor(as.numeric(DF$macro))))
collectGarbage()
nMeta = exprSize$nGenes
nSamples = exprSize$nSamples

save(multiExpr, Traits, nMeta, nSamples, setLabels, shortLabels, exprSize, file = "Consensus-dataInput.RData")
rm(list=ls())
####---------- Step 2: network construction and module detection ----------####
#### ! The one-step automatic function cannot do weighted network analysis
#### We use step-by-step method here for consensus network analysis
## Setup the R session
enableWGCNAThreads()
set.seed(2)
lnames = load(file = "Consensus-dataInput.RData")
# The variable lnames contains the names of loaded variables
lnames
nSets = checkSets(multiExpr)$nSets

## Network construction and module detection
# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by = 1), seq(12,30,by = 2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for(set in 1:nSets){
  weights = matrix(rep(1/Traits[[set]]$data$v, nMeta), nSamples[set], nMeta)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector = powers,
                                                     weights = weights,networkType = "signed",
                                                     verbose = 2)[[2]])
}

collectGarbage()
# plot the results
colors = c("black","red","blue","green")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8,6)
pdf(file = "Figures/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

### Network construction and consensus module detection

## Calculation of network adjacencies
softPower = c(16,24,30,24)

# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0,dim = c(nSets, nMeta, nMeta))

# Calculate adjacencies in each individual dataset
for(set in 1:nSets){
  weights = 1/Traits[[set]]$data$v
  adjacencies[set,,] = ((cor(multiExpr[[set]]$data,use="p", weights.x = weights) + 1)/2)^softPower[set]
}

## Calculation of Topological overlap
#Initialize an appropriate array to hold the TOMs
TOM = array(0, dim=c(nSets, nMeta, nMeta))
for(set in 1:nSets){
  TOM[set,,] = TOMsimilarity(adjacencies[set,,])
}

## Scaling of Topological Overlap matrices to make them comparable across sets
# Define the reference percentile
scaleP = 0.95
#set RNG seed for reproducibility of sampling
set.seed(12345)
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nMeta * (nMeta - 1)/2, size = nSamples)
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1,nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
for(set in 1:nSets){
  TOMScalingSamples[[set]] = as.dist(TOM[set,,])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8)
  # Scale the other visits' TOM
  if(set > 1){
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
    TOM[set,,] = TOM[set,,] ^ scalePowers[set]
  }
}

scaledTOMSamples = list()
for(set in 1:nSets){
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
}

sizeGrWindow(6,8)
pdf(file = "Figures/TOMScaling-QQPlot.pdf", wi = 8, he = 6)
par(mfrow=c(1,3))
#qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = T,
                    cex = 0.6, xlab = paste("TOM in", setLabels[1]),
                    ylab = paste("TOM in", setLabels[2]), main = "Q-Q plot of TOM",
                    pch = 20)
#qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = F)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a = 0, b = 1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black","red"))

qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[3]], plot.it = T,
                    cex = 0.6, xlab = paste("TOM in", setLabels[1]),
                    ylab = paste("TOM in", setLabels[3]), main = "Q-Q plot of TOM",
                    pch = 20)
#qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[3]], plot.it = F)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a = 0, b = 1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black","red"))

qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[4]], plot.it = T,
                    cex = 0.6, xlab = paste("TOM in", setLabels[1]),
                    ylab = paste("TOM in", setLabels[4]), main = "Q-Q plot of TOM",
                    pch = 20)
#qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[4]], plot.it = F)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a = 0, b = 1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black","red"))
dev.off()

## Calculation of consensus Topological Overlap
consensusTOM = pmin(TOM[1,,],TOM[2,,],TOM[3,,], TOM[4,,])
## Clustering and module identification
# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")
minModuleSize = 15
# Module identifcation using dynamic tree cut
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,pamRespectsDendro = F)
unmergedColors = labels2colors(unmergedLabels)
table(unmergedLabels)
sizeGrWindow(8,6)
pdf(file = "Figures/ConsensusDendrogram-unmerged.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree,unmergedColors,"Dynamic Tree Cut",
                    dendroLabels= F, hang = 0.03, addGuide = T, guideHang = 0.05)
dev.off()

## Merging of modules whose expression profiles are very similar
# Calculate module eigen metabolites
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigen metabolites
consMEDiss = consensusMEDissimilarity(unmergedMEs)
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average")
# plot the results
sizeGrWindow(7,6)
par(mfrow=c(1,1))
plot(consMETree,main = "Consensus clustering of consensus module eigen metabolites",
     xlab = "", sub = "")
abline(h = 0.25,col = "red")

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors
moduleColors = labels2colors(moduleLabels)
table(moduleLabels)
table(moduleColors)

#save a copy:
modMemb <- as.data.frame(cbind(moduleColors, moduleLabels))
modMemb$met_ID <- substring(names(multiExpr[[1]]$data),4)
modMemb <- inner_join(modMemb,anno,by="met_ID", y.all = F, x.all=T)
mm.tab <- table(modMemb$Lipid.class, modMemb$moduleColors)
View(unclass(mm.tab))
write.csv(modMemb, "Result/network/Module membership_cons.csv", row.names=F)
# Eigen metabolites of the new merged modules:
consMEs = merge$newMEs
# plot the dendrogram:
sizeGrWindow(9,6)
plotDendroAndColors(consTree,cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
save(consMEs, moduleColors, moduleLabels, consTree, file = "Consensus-NetworkConstruction-man.RData")


####--------------- Step 3: Relating consensu modules to set-specific modules ---------------####
lnames = load(file = "Consensus-dataInput.RData")
lnames
lnames = load(file = "Consensus-NetworkConstruction-man.RData")
lnames
lnames = load(file = "Visit-4-networkConstruction-man.RData")
lnames
Visit0Labels = moduleLabels
Visit0Colors = moduleColors
Visit0Tree = geneTree
Visit0MEs = orderMEs(MEs,greyName = "ME0")

lnames = load(file = "Consensus-NetworkConstruction-man.RData")
lnames

# Isolate the module labels in the order they appear in ordered module eigengenes
Visit0Modules = substring(names(Visit0MEs),3)
consModulelabels = substring(names(consMEs[[1]]$data),3)

# Convert the numeric module labels to color labels
consModules = labels2colors(as.numeric(consModulelabels))
label_color <- data.frame(label = as.integer(consModulelabels), color = consModules)
# Number of Visit 0 and consensus modules
nVisit0Mods = length(Visit0Modules)
nConsMods = length(consModules)

# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0,nVisit0Mods,nConsMods)
countTbl = matrix(0, nVisit0Mods, nConsMods)

# Execute all pairwise comparisons
for(vmod in 1:nVisit0Mods){
  for(cmod in 1:nConsMods){
    VisMembers = (Visit0Colors == Visit0Modules[vmod])
    consMembers = (moduleColors == consModules[cmod])
    pTable[vmod, cmod] = -log10(fisher.test(VisMembers, consMembers, alternative = "greater")$p.value)
    countTbl[vmod, cmod] = sum(Visit0Colors == Visit0Modules[vmod] & moduleColors == consModules[cmod])
  }
}

## plot the result
# Truncate p-values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3 * max(pTable[is.finite(pTable)])
pTable[pTable > 50] <- 50
# Margininal counts (really module sizes)
Visit0ModTotals = apply(countTbl,1,sum)
consModTotals = apply(countTbl,2,sum)
# Actual plotting
sizeGrWindow(3.03,1.99)
pdf(file ="~/Collaboration_fetal/lipidomics/version 6.0/Cons-Visit4.pdf", 
    wi = 3.03, he = 1.99)
par(mfrow=c(1,1))
par(cex=0.4)
par(mar=c(6.5,7.5,1.5,0.5))
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", Visit0Modules),
               colorLabels = T,
               xSymbols = paste("Cons ", consModules,": ", consModTotals, sep = ""),
               ySymbols = paste("Visit 4 ", Visit0Modules, ": ", Visit0ModTotals, sep = ""),
               textMatrix = countTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of Visit 4-specific and consensus modules",
               cex.text = 1.2, cex.lab = 0.9, setStdMargins = F)
dev.off()


####--------------- Step 4: Relating consensus modules to sample information ---------------####
# Load the data saved in the first part
lnames = load("Consensus-dataInput.RData")
lnames = load("Consensus-NetworkConstruction-man.RData")
lnames
exprSize = checkSets(multiExpr)
nSets = exprSize$nSets

# Set uo variables to contain the module-trait correlations
moduleTraitCor = vector(mode = "list", length = nSets)
moduleTraitPvalues = vector(mode = "list", length = nSets)

# Calculate the correlations
for(set in 1:nSets){
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data[,-c(1:3)], 
                              weights.y = 1/Traits[[set]]$data$v, use = "p")
  moduleTraitPvalues[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
}

# Convert numerical labels to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data),3)))
MEColorName = paste("ME", MEColors)

# plot the module-trait relationship table for set number 1
for(set in 1:nSets){
  pdf(file=paste("Figures/ModuleTraitRelationships-",shortLabels[set], ".pdf",sep = ""), wi = 10, he = 7)
  textMatrix = paste(signif(moduleTraitCor[[set]],2), "\n(",
                   signif(moduleTraitPvalues[[set]],1),")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor[[set]])
  par(mar=c(8,8.8,3,2.2))
  labeledHeatmap(Matrix = moduleTraitCor[[set]],
                 xLabels = names(Traits[[set]]$data[,-c(1:3)]),
                 yLabels = MEColorName,
                 ySymbols = MEColorName,
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module--trait relationship in", setLabels[set]))
  dev.off()

}
# for each module-trait pair, we take 
# the correlation that has the lower absolute valaue in the four sets if they have the same sign
# and zero if one have opposite signs.
# Initialize matrices to hold the consensus correlation and p-value

consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalues = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))

# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 & moduleTraitCor[[4]] < 0
consensusPvalues[negative] = pmax(moduleTraitPvalues[[1]][negative], moduleTraitPvalues[[2]][negative], 
                                  moduleTraitPvalues[[3]][negative], moduleTraitPvalues[[4]][negative])
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                              moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative])

# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 & moduleTraitCor[[4]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                              moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive])
consensusPvalues[positive] = pmax(moduleTraitPvalues[[1]][positive], moduleTraitPvalues[[2]][positive],
                                 moduleTraitPvalues[[1]][positive], moduleTraitPvalues[[1]][positive])
textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalues, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Figures/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data[,-c(1:3)]),
               yLabels = MEColorName,
               ySymbols = MEColorName,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels[1:3], collapse = ", "), "and ",setLabels[4]))
dev.off()

## Exporting results of the network analysis
probes = data.frame(met_ID = substring(names(multiExpr[[1]]$data),4))
probes = inner_join(probes, anno[,c("met_ID","Annotation")])
probes$Ann. = substring(probes$Annotation,1,3)

## recalculate the module eigengenes in the "alphabetic" order and 
## calculate the gene significance and module memberships in each data set
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels)
kME = list()
for(set in 1:nSets){
  weights = 1/Traits[[set]]$data$v
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data, weights.x = weights)
}
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z + kME[[3]]$Z + kME[[4]]$Z)/sqrt(4)
kME.metaP = 2 * pnorm(abs(kME.metaZ), lower.tail = F)

kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[3]]$cor, kME[[4]]$cor,
               kME[[1]]$p, kME[[2]]$p, kME[[3]]$p, kME[[4]]$p,
               kME.metaZ, kME.metaP)
MEnames = colnames(consMEs.unord[[1]]$data)
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nMeta, 10*nMEs)
rownames(kMEmat) = probes$met_ID
colnames(kMEmat) = spaste(c("kme.set1","kme.set2","kme.set3","kme.set4",
                            "p.kme.set1","p.kme.set2","p.kme.set3","p.kme.set4",
                            "Z.kME.meta", "p.kME.meta"), rep(MEnames,rep(10,nMEs)))
info = data.frame(probes,ModuleLabel = moduleLabels, ModuleColor = labels2colors(moduleLabels), kMEmat)
write.csv(info, file = "Result/module-info.csv", row.names = F)

#### ------------------- Step 5: Compare eigengene networks in different visits --------------- ####
anthro = vector(mode = "list", length = nSets)
for(set in 1:nSets){
  anthro[[set]] = list(data=as.data.frame(Traits[[set]]$data$LGA)) #, Traits[[set]]$data$bwz,Traits[[set]]$data$Length,
                                          #Traits[[set]]$data$SSF, Traits[[set]]$data$head_circ,Traits[[set]]$data$ab_circ))
  names(anthro[[set]]$data) = c("LGA")#,"bwz","len","SSF","head_circ","ab_circ")
}

# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)

# add the trait to the eigengenes and order them by consensus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC,anthro))
MET_woanthro = consensusOrderMEs(consMEsC)
# perform the differential analysis with plotEigengeneNetworks
sizeGrWindow(16,16.5)
pdf(file = "Figures/network/EigengeneNetworks.pdf", width= 16, height = 16.5);
par(cex = 0.8)
plotEigengeneNetworks(MET_woanthro, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90, excludeGrey = F,
                      greyLabel = "black")

dev.off()

# Set up the function of fitting anthro~ME relationship
moduleTraitBeta <- moduleTraitPvalues <- moduleTraitPFDR <- vector(mode = "list", length = nSets)
lme_ME <- function(Cov, idx, x){
  nDF <- cbind(Cov[,c(1:10,idx)], x)
  nDF <- na.omit(nDF)
  names(nDF)[11] = "anthro"
  model = lmer(anthro ~ x + age + race + education + preBMI + GA_blood + factor(Gender) + nulli +( 1 | pairID), 
              data = nDF,weights = 1/v)
  est <- summary(model)$coefficients[2,]
  return(est)
}

lme_ME_binary <- function(Cov, idx, x){
  nDF <- cbind(Cov[,c(1:10,idx)], x)
  nDF <- na.omit(nDF)
  names(nDF)[11] = "anthro"
  model = glmer(anthro ~ x + age + race + education + preBMI + GA_blood + factor(Gender) + nulli +( 1 | pairID), 
               data = nDF,weights = 1/v, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10, family = binomial)
  est <- summary(model)$coefficients[2,]
  return(est)
}

ind_test <- NULL
idx = 20
for(i in 1:nSets){
  Cov <- Traits[[i]]$data
  #Cov$LGA = Cov$LGA/1000
  MetMEs <- consMEs[[i]]$data
  res <- as.data.frame(t(sapply(MetMEs, lme_ME_binary, Cov = Cov, idx = idx)))[,c(1,2,4)]
  tab.res <- data.frame(label = as.numeric(substring(rownames(res),3)), 
                        beta_se = paste0(round(res[,1],4)," (", round(res[,2],4), ")"),
                        beta = res[,1],
                        se = res[,2],
                        pvalue = res[,3])
  tab.res$pFDR <- p.adjust(tab.res$pvalue, "fdr")
  tab.res <- inner_join(tab.res, label_color, by = "label")
  tab.res$Visit = ifelse(i == 4, i, i-1)
  ind_test <- rbind(ind_test, tab.res[,-2])
  
  moduleTraitBeta[[i]] <- tab.res[,c("label", "color", "beta_se")]
  moduleTraitPvalues[[i]] <- tab.res[,c("label", "color", "pvalue")]
  moduleTraitPFDR[[i]] <- tab.res[,c("label", "color", "pFDR")]
}



write.csv(ind_test, file = paste0('Result/network/',names(Cov)[idx],"_cons.csv"), row.names = F)

 ###-----------------Multiple regression ----------------
for(idx in 11:16){
  multi_test <- NULL
  for(i in 1:nSets){
    Cov <- Traits[[i]]$data
    Cov$bw = Cov$bw/1000
    MetMEs <- consMEs[[i]]$data
    nME = ncol(MetMEs)
    nDF <- cbind(Cov[,c(1:10,idx)], MetMEs)
    nDF <- na.omit(nDF)
    names(nDF)[11] = "anthro"
    model = lmer(as.formula(paste("anthro ~", paste(names(MetMEs), collapse = "+"),
                                  "+ age + race + education + preBMI + GA_blood + factor(Gender) + nulli +( 1 | pairID)")), 
                 data = nDF,weights = 1/v)
    
    res <- summary(model)$coefficients[1+(1:nME),c(1,2,5)]
    
    tab.res <- data.frame(label = as.numeric(substring(rownames(res),3)), 
                          beta = res[,1],
                          se = res[,2],
                          pvalue = res[,3])
    tab.res <- inner_join(tab.res, label_color, by = "label")
    tab.res$Visit = ifelse(i == 4, i, i-1)
    multi_test <- rbind(multi_test, tab.res[,-2])
  }
  
  write.csv(multi_test, file = paste0('Result/network/',names(Cov)[idx],"_cons_multipleRegression.csv"), row.names = F)
  
}

# Convert numerical labels to colors for labeling of modules in the plot

textMatrix = paste(signif(moduleTraitBeta[[set]]$beta,2),"\n(",
                   signif(moduleTraitPvalues[[set]]$pvalue,1), ")", sep="")
for(set in 2:nSets){
  cbind(textMatrix, paste(signif(moduleTraitBeta[[set]]$beta,2),"\n(",
                          signif(moduleTraitPvalues[[set]]$pvalue,1), ")", sep=""))
}

###5.Intramodular analysis: Identify top metabolites for each module --------------------------------------
i = 0
for(visit in c(0,1,2,4)){
  i = i+1
  #Module color and sample size
  DFmet <- multiExpr[[i]]$data
  modNames = label_color[label_color[,1] == as.integer(substring(names(consMEs[[i]]$data), 3)),2]
  nSample = nrow(DFmet)

  #Module membership (MM): correlation btw metabolites and module eigengene
  MM   <- as.data.frame(cor(DFmet, consMEs[[i]]$data, use = "p"))
  p.MM <- as.data.frame(corPvalueStudent(as.matrix(MM), nSample))
  names(MM)   <- paste("MM", modNames, sep="")
  names(p.MM) <- paste("p.MM", modNames, sep="")
  MM <- cbind(MM, p.MM)  
  MM$lipid_identifier   <- rownames(MM)
  MM$met_ID <- substring(MM$lipid_identifier,4)
  #Gene Significance(GS): the abs. correlation btw each metabolite and outcome
  GS <- read.csv("tmp/known/bw1.csv", header=T, sep=",")
  GS <- GS[,c("met_ID", "Lipid.class","Annotation", "beta", "pval","pFDR")]
  GS <- inner_join(GS,modMemb[,c("met_ID","moduleColors","moduleLabels")],by = "met_ID")
  GS$anno_sig <- GS$Annotation
  GS[GS$pFDR>=0.05, "anno_sig"] <- ""     #only show label for p < 0.05/327 (Bonferroni)
  GS$metID_sig <- GS$met_ID
  GS[GS$pFDR>=0.05, "metID_sig"] <- ""  #only show label for p < 0.05/327 (Bonferroni)

  #Merge GS and MM - save a copy
  MM.GS <- merge(MM, GS, "met_ID")
  write.csv(MM.GS, file=paste("Result/network/bw_intramodular_MM vs. GS_cons", visit, ".csv", sep=""))

#Scatterplot of Gene Significance vs. Module Membership in the "brown" or "turquoise" module (highest corr with GDM)

for (module in modNames) {
  row <- (MM.GS$moduleColors==module)
  png(paste("Intramodular_MM vs GS_signed_cons", visit, "_", module, ".png", sep=""), width=4800, height=4800, res=500)
  #sizeGrWindow(4800,4800)
  verboseScatterplot(MM.GS[row, paste("MM",module,sep="")],
                     MM.GS[row, "beta"],
                     abline = TRUE,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Beta for Birthweight Z-score" ,
                     main = paste("MM vs GS, visit", visit, ", ", module),
                     cex.main = 1, cex.lab = .8, cex.axis = .8, col = module, pch=16)
  text(MM.GS[row, paste("MM",module,sep="")], 
       MM.GS[row, "beta"], 
       labels = MM.GS[row, "anno_sig"], cex=0.8, pos=2)
  dev.off()
}
}
# Meaningful if GS and MM are highly correlated
# after combining with annotation file, can 
# Order the metabolites first by module color, then by geneTraitSignificance
# to look out central players in each module

#### trajectory analysis of modules:
lnames = load("Consensus-dataInput.RData")
lnames = load("Consensus-NetworkConstruction-man.RData")
lnames
nME = ncol(consMEs[[1]]$data)
label_color <- unique(cbind(moduleLabels,moduleColors))

for(j in 1:nME){
  consME1 <- data.frame(ID = row.names(consMEs[[1]]$data), ME1 = consMEs[[1]]$data[,j])
  consME2 <- data.frame(ID = row.names(consMEs[[2]]$data), ME2 = consMEs[[2]]$data[,j])
  consME3 <- data.frame(ID = row.names(consMEs[[3]]$data), ME3 = consMEs[[3]]$data[,j])
  consME4 <- data.frame(ID = row.names(consMEs[[4]]$data), ME4 = consMEs[[4]]$data[,j])
  consME <- merge(merge(merge(consME1, consME2, by = "ID", all = T), consME3, by = "ID", all = T),
                  consME4, by = "ID", all = T)
  traits1 <- data.frame(ID = Traits[[1]]$data$ID, wt = 1/Traits[[1]]$data$v, 
                        T1 = Traits[[1]]$data$GA_blood, TCOV1 = Traits[[1]]$data$diet)
  traits2 <- data.frame(ID = Traits[[2]]$data$ID, wt = 1/Traits[[2]]$data$v,
                        T2 = Traits[[2]]$data$GA_blood, TCOV2 = Traits[[2]]$data$diet)
  traits3 <- data.frame(ID = Traits[[3]]$data$ID,
                        T3 = Traits[[3]]$data$GA_blood, TCOV3 = Traits[[3]]$data$diet)
  traits4 <- data.frame(ID = Traits[[4]]$data$ID, 
                        T4 = Traits[[4]]$data$GA_blood, TCOV4 = Traits[[4]]$data$diet)
  traits <- merge(merge(merge(traits1, traits2, by = c("ID","wt"), all = T), traits3, by = "ID", all = T),
                            traits4, by = "ID", all = T)
  data <- merge(consME,traits,by = "ID", all = T)
  data <- data[,c(1,6,2:5,7,9,11,13,8,10,12,14)]
  write.csv(data, file = paste0('Data/sas_data/cons/',names(consMEs[[1]]$data)[j],".csv"), row.names = F)
}

#### trajectory vs. outcomes ####
library(tidyverse)
library(forcats)
library(lme4)
library(lmerTest)
library(car)
anthro <- read.csv("traj_anthro.csv")
anthro$education <- factor(anthro$education > 3, ordered = T)
anthro$race <- factor(anthro$race == 3)
anthro$nulli <- factor(anthro$nulli)
anthro$alcohol <- factor(anthro$alcohol)
anthro$Gender <- factor(anthro$Gender)
#anthros = c("bw","bwz","Length","head_circ","ab_circ","SSF") #"LGA","SGA","macro","LBW"
anthros="LGA"
anthros_idx <- which(names(anthro) %in% anthros)
anthros <- names(anthro)[anthros_idx]

res <- NULL
for(j in 1:6){
  idx = anthros_idx[j]
  pval <- rep(NA, 7)
  for(i in 0:6){
    DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                     age = anthro$age, race = anthro$race, education = anthro$education, 
                     preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                     outcome = anthro[,idx])
    
    grouping <- read.csv(paste0("Data/sas_data/cons/res/out/outME",i,".csv"))
    DF <- merge(DF, grouping[,c("ID","GROUP")],by="ID")
    DF$GROUP <- fct_infreq(as.factor(DF$GROUP))
  
    DF <- na.omit(DF)
    model = glmer(outcome ~ factor(GROUP) + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
                  data = DF, family = binomial)
    #model = glmer(outcome ~ factor(GROUP) + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
    #             data = DF, weights = wt, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
    #model1 = lmer(anthro ~ factor(GROUP) + age + race + education + preBMI + GA_blood + Gender + nulli + alcohol + post + (1 | pairID), 
    #            data = DF, weights = 1/v)
    pval[i+1] <- Anova(model,type=2)[[3]][1]
  }
  traj_res <- data.frame(module=paste0("ME",0:6), pval = pval, pBonf = p.adjust(pval,"bonferroni"), pFDR = p.adjust(pval,"fdr"))
  traj_res$moduleLabels <- substring(traj_res$module,3)
  traj_res <- merge(traj_res,label_color,by="moduleLabels")
  traj_res <- traj_res[,c(1:2,6,3:5)]
  traj_res$outcome <- anthros[j]
  
  res <- rbind(res,traj_res)
}
write.csv(res, file = "Result/network/Module traj_test.csv", row.names = F)

#### Individual trajectory analysis
# Length vs. grey and red
# be vs. grey

##### -------------- step 1: plot out the trajectories --------------
os <- read.csv("Data/sas_data/cons/res/os/osME0.csv")
op <- read.csv("Data/sas_data/cons/res/op/opME0.csv")

plot.data <- data.frame(Time = rep(op$T,4),
                        Prediction = unlist(op[,6:9]),
                        Beta0 = rep(os$BETA0,each = 4),
                        Upper = unlist(op[,c(11,13,15,17)]),
                        Lower = unlist(op[,c(10,12,14,16)]),
                        Pattern = rep(paste(1:4),each = 4),
                        Perc = rep(os$PI,each = 4))
plot.data$label <- paste("Class",plot.data$Pattern, round(plot.data$Perc,2),"%")
ggplot(plot.data,aes(x = Time, col = label))+
  geom_pointrange(aes(y = Prediction, ymin = Lower, ymax= Upper), linetype = 2)+
  geom_line(aes(y = Prediction)) + 
  theme_classic()+
  labs(x = "Gestational Weeks", y = "Predicted Scores on Grey Module")+
  scale_color_discrete(name = "Trajectory Class")+
  theme(legend.position = "bottom")

os <- read.csv("Data/sas_data/cons/res/os/osME6.csv")
op <- read.csv("Data/sas_data/cons/res/op/opME6.csv")

plot.data <- data.frame(Time = rep(op$T,3),
                        Prediction = unlist(op[,5:7]),
                        Upper = unlist(op[,c(9,11,13)]),
                        Lower = unlist(op[,c(8,10,12)]),
                        Pattern = rep(paste(1:3),each = 4),
                        Perc = rep(os$PI,each = 4))
plot.data$label <- paste("Class",plot.data$Pattern, round(plot.data$Perc,2),"%")
ggplot(plot.data,aes(x = Time, col = label))+
  geom_pointrange(aes(y = Prediction, ymin = Lower, ymax= Upper), linetype = 2)+
  geom_line(aes(y = Prediction)) + 
  theme_classic()+
  labs(x = "Gestational Weeks", y = "Predicted Scores on Red Module")+
  scale_color_discrete(name = "Trajectory Class")+
  theme(legend.position = "bottom")
library(emmeans)
DF <- data.frame(ID = anthro$ID, pairID = anthro$pairID, wt = anthro$wt,
                age = anthro$age, race = anthro$race, education = anthro$education, 
                preBMI = anthro$preBMI, nulli = anthro$nulli, Gender = anthro$Gender,
                outcome = anthro$bw)
grouping <- read.csv("Data/sas_data/cons/res/out/outME0.csv")
DF <- merge(DF, grouping[,c("ID","GROUP")],by="ID")
DF$GROUP <- fct_infreq(as.factor(DF$GROUP))

DF <- na.omit(DF)
model = lmer(outcome ~ GROUP + age + race + education + preBMI + Gender + nulli+ (1 | pairID), 
             data = DF, weights = wt)
model.anova <- Anova(model,type = "II")
multi <- emmeans(model, list(pairwise ~ GROUP), adjust = "tukey")
save(multi,model,file = "Result/network/multi_grey_bw.Rdata")

# read in the trajectory analysis results
indi_traj <- read.csv("Result/trajectory analysis_continuous/table3_length.csv")
cons_mem <- read.csv("Result/Module membership_cons.csv")
indi_traj <- merge(indi_traj, cons_mem[,c("met_ID","moduleColors")],by = "met_ID", all.x = T)
sig_traj <- indi_traj[indi_traj$moduleColors=="red",c("met_ID","Annotation","Lipid.class","pFDR")]
write.csv(sig_traj,'Result/network/sig_cons_red_length.csv',row.names = F)

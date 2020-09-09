library(WGCNA)
library(tidyverse)
setwd("/home/scott/Documents/Uni/Research/Projects/Blueberry_Data/WGCNA_Data")
options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)
allowWGCNAThreads(6)
WGCNAnThreads()

# Read in the blueberry data set
# Each row is a gene and each column is a library.
BlueberryData = read.csv("SingleHaplotype_Blueberry_TPM.csv", header=TRUE, row.names='Gene_Name')


the_genes = rownames(BlueberryData)
the_libraries = colnames(BlueberryData)
n_the_libraries = ncol(BlueberryData)
# Subset data based on whether a single entry is greater than 10, take that row
# 10 for testing purposes, 0.1 for general guidelines
subset_val = 0.1
BlueberryData = filter_all(as_tibble(BlueberryData), any_vars(. > subset_val))
BlueberryData = t(BlueberryData)

# Decide what soft-threshold to use
soft_threshold_graph = function(){
        # Choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        # Call the network topology analysis function
        sft = pickSoftThreshold(BlueberryData, powerVector=powers, verbose = 5, blockSize=15000)
        # Plot the results:
        sizeGrWindow(9, 5)
        par(mfrow = c(1,2));
        cex1 = 0.9;
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red");
        # this line corresponds to using an R^2 cut-off of h
        abline(h=0.90,col="red")
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

# Choose a soft power of 6 because that was best
softPower = 6
adjacency = adjacency(BlueberryData, power=softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
#rm(adjacency)
rm(TOM)

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
gtree_name = 'RawGeneTree_'
subset_val = toString(subset_val)
#sizeGrWindow(12,9)
ext = '.pdf'
pdf(file=paste(gtree_name, subset_val, ext,  sep=''), width=12, height=10)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
#0.04, look up guidehang
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
			    deepSplit = 2, pamRespectsDendro = FALSE,
			    minClusterSize = minModuleSize)
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
gtree_name = 'ColorGeneTree_'
subset_val = toString(subset_val)
ext = '.pdf'
pdf(paste(gtree_name, subset_val, ext,  sep=''))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
		    dendroLabels = FALSE, hang = 0.03,
		    addGuide = TRUE, guideHang = 0.05,
		    main = "Gene dendrogram and module colors")
dev.off()
table(dynamicColors)

# Calculate the eigengenes
# Merging of modules whose expression profiles are similar
MEList = moduleEigengenes(BlueberryData, colors=dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
#Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')
gtree_name = 'MergeColorGeneTree_'
pdf(paste(gtree_name, subset_val, ext,  sep=''))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


# Plot the cut line into the dendrogram
# Corresponding to correlation of 0.75
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(BlueberryData, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
dev.off()

gtree_name = 'Merged_and_Raw_ColorGeneTree_'
pdf(paste(gtree_name, subset_val, ext,  sep=''))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(100))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
rm(mergedMEs)
# Save module colors and labels for use in subsequent parts
save(MEs, MEList, moduleLabels, moduleColors, geneTree, BlueberryData, file = 'Blueberry_GeneNetwork.RData')
load('Blueberry_GeneNetwork.RData')
names(BlueberryData)

# --------------------------
# Tutorial 3
modNames = substring(names(MEs), 3)
module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module



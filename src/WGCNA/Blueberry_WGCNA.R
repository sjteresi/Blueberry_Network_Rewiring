# __author__ = Scott Teresi
library(WGCNA)
library(tidyverse)
library(docstring)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)
allowWGCNAThreads(6)
WGCNAnThreads()


load_all_hap_hpcc = function() {
  # Read in the blueberry data set
  # Each row is a gene and each column is a library.
  setwd("/mnt/research/edgerpat_lab/Scotty/Blueberry_Network_Rewiring/results/FPKM_TPM")
  AllHap_BlueberryData = read.csv("Blueberry_TPM_All_Haplotype.tsv", header = TRUE, sep =
                                    '\t')
  return(AllHap_BlueberryData)
}

soft_threshold_graph = function(data_matrix) {
  # Generate a graph displaying the scale-free topology fit index as
  # a function of soft-thresholding power. Only needs to be run once for a given dataset,
  # Pick the first value that is over the red cut line. Lower values are better.
  # Args:
  # 	data_matrix (samples, genes): TPM matrix
  # Returns:
  #	None.
  #	Generates graph. Saves to file
  
  # Decide what soft-threshold to use
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(
    data_matrix,
    powerVector = powers,
    verbose = 5,
    blockSize = 15000
  )
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  
  cex1 = 0.9
  
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  print('NOTE: Soft thresholding graph')
  pdf(file = 'Scale_Free_Topology_Fit.pdf',
      width = 8,
      height = 5)
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  dev.off()
  print('NOTE: Mean connectivity graph')
  pdf(file = 'Mean_Connectivity.pdf',
      width = 8,
      height = 5)
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
}

calc_diss_TOM = function(adjacency) {
  # Transform adjacency matrix into Topological Overlap Matrix and calculate dissimilarity
  # Args:
  # adjacency ():
  # Returns:
  # dissTOM ():
  print("NOTE: Calculating TOM from adjacency")
  TOM = TOMsimilarity(adjacency)
  write.table(
    TOM,
    'ModuleMembership_pvalue.tsv',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
  )
  # NB had some correspondence with Peter Langfelder and came to the conclusion that
  # there was an error in the matrix math libraries (was having an issue between
  # HPCC and local machine). Matrix math library for HPCC now specified in job submission script.
  # This is why it is always good to read the documentation well because my TOM and dissTOM
  # values were unexpected and it wasn't breaking the code or anything.

  # print("checking peter's suggestion")
  # test_val = max(TOM - 1)
  # print(test_val)
  dissTOM = 1 - TOM
  return(dissTOM)
}

cluster_using_TOM = function(dissTOM) {
  # Now we cluster using TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")

  # Plot the resulting clustering tree (dendrogram)
  gtree_name = 'RawGeneTree'
  ext = '.pdf'
  print('NOTE: Cluster using TOM')
  pdf(
    file = paste(gtree_name, ext,  sep = ''),
    width = 8,
    height = 6
  )

  plot(
    geneTree,
    xlab = "",
    sub = "",
    main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE,
    hang = 0.04
  )
   # MAGIC hang value

  dev.off()
  return(geneTree)
}


trim_and_color = function(geneTree, dissTOM) {
  # We like large modules, so we set the minimum module size relatively high:
  # MAGIC minModuleSize
  minModuleSize = 30
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize
  )
  
  
  # Convert numeric lables into colors
  dynamicColors <<- labels2colors(dynamicMods)

  # Plot the dendrogram and colors underneath
  gtree_name = 'ColorGeneTree'
  ext = '.pdf'
  print('NOTE: Trim and color ')
  pdf(paste(gtree_name, ext,  sep = ''),
      height = 6,
      width = 8)
  plotDendroAndColors(
    geneTree,
    dynamicColors,
    "Dynamic Tree Cut",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
  dev.off()
  #table(dynamicColors)
  
}

merging_similar_expression = function(working_data, dynamicColors, geneTree) {
  # Calculate the eigengenes
  # Merging of modules whose expression profiles are similar
  MEList = moduleEigengenes(working_data, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1 - cor(MEs)

  #Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = 'average')
  gtree_name = 'MergeColorGeneTree'
  ext = '.pdf'
  print('NOTE: Merging similar expression')
  pdf(paste(gtree_name, ext,  sep = ''),
      width = 9,
      height = 6)
  plot(METree,
       main = "Clustering of module eigengenes",
       xlab = "",
       sub = "")
  
  # Plot the cut line into the dendrogram
  # Corresponding to correlation of 0.75
  # NOTE MAGIC number
  MEDissThres = 0.25
  abline(h = MEDissThres, col = "red")

  # Call an automatic merging function
  merge = mergeCloseModules(working_data,
                            dynamicColors,
                            cutHeight = MEDissThres,
                            verbose = 3)

  # The merged module colors
  mergedColors = merge$colors

  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  dev.off()
  
  gtree_name = 'Merged_and_Raw_ColorGeneTree'
  #jpeg(paste(gtree_name, ext,  sep=''), width=800, height=600)
  print('NOTE: Merging similar expression new graph')
  pdf(paste(gtree_name, ext,  sep = ''),
      width = 8,
      height = 6)
  plotDendroAndColors(
    geneTree,
    cbind(dynamicColors, mergedColors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  # Rename to moduleColors
  moduleColors <<- mergedColors

  # Construct numerical labels corresponding to the colors
  colorOrder <<- c("grey", standardColors(90))
  moduleLabels <<- match(moduleColors, colorOrder) - 1
  MEs <<- mergedMEs

  # Save module colors and labels for use in subsequent parts
  # NOTE, this corresponds to a checkpoint in the tutorial
  save(MEs, moduleLabels, moduleColors, geneTree, file = "Chapter2.RData")
  
}

check_bad_genes = function(expression_dataframe, test_type) {
  gsg = goodSamplesGenes(expression_dataframe, verbose = 3)
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes) > 0)
      sink(
        paste(test_type, 'removed_genes.txt', sep = '_'),
        append = FALSE,
        split = FALSE
      )
    printFlush(paste("Removing genes:\n", paste(
      names(expression_dataframe)[!gsg$goodGenes], collapse = "\n"
    )))
    sink()
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples.txt", paste(
        rownames(expression_dataframe)[!gsg$goodSamples], collapse = ", "
      )))
    
    # Remove the offending genes and samples from the data:
    expression_dataframe = expression_dataframe[gsg$goodSamples, gsg$goodGenes]
  }
  return(expression_dataframe)
  
}

names_and_colors = function(working_data) {
  ngenes = ncol(working_data)
  nsamples = nrow(working_data)
  MEsO = moduleEigengenes(working_data, moduleColors)$eigengenes
  MEs = orderMEs(MEsO)
  # ^ Module eigenvectors that are ordered

  modNames = substring(names(MEs), 3)
  geneModuleMembership <<-
    as.data.frame(cor(working_data, MEs, use = 'p'))
  # ^ computes the covariance, p is just a method for handling missing values, not a p-value thing.
  
  MMPvalue <<-
    as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nsamples))
  names(geneModuleMembership) = paste("MM_", modNames, sep = '')
  names(MMPvalue) = paste('p_MM_', modNames, sep = '')
  
  mmpvalue_dframe = tibble::rownames_to_column(MMPvalue, 'Gene_Names')
  write.table(
    mmpvalue_dframe,
    'ModuleMembership_pvalue.tsv',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
  )
  
 
  # NOTE: 
  # How to get thegenes for a particular module
  #module = 'thistle'
  #column = match(module, modNames)
  #moduleGenes = moduleColors == module
  #names(working_data)
}

process = function(load_data_function, test_type, softpower) {
  # Process data, master control function
  #
  #  Args:
  #      load_data_function (function): A function that returns a TPM matrix
  #      softpower (int): An integer representing the soft threshold value to pick to satisfy
  #           scale-free topology, investigate the soft_threshold_graph for more.
  #           This only needs to be solved once.
  
  
  
  BlueberryData = load_data_function  # global variable for inspection

  # Set output directory, MAGIC HPCC path
  setwd("/mnt/research/edgerpat_lab/Scotty/Blueberry_Network_Rewiring/results/WGCNA")
  transposed_BlueberryData = as.data.frame(t(subset(BlueberryData, select = -c(Gene_Name))))
  names(transposed_BlueberryData) = BlueberryData$Gene_Name
  rownames(transposed_BlueberryData) = names(subset(BlueberryData, select = -c(Gene_Name)))
  
  # data should be in shape (samples, genes) before use
  working_data <<-
    check_bad_genes(transposed_BlueberryData, test_type)
  rm(BlueberryData, transposed_BlueberryData)
   
  # NOTE run this if you have not previously examined the network topology
  # and picked a soft thresholding power, command may take some time, and a lot of RAM
  soft_threshold_graph(working_data)
  
  print('NOTE: Calculating adjacency values, may take time')
  adjacency = adjacency(working_data, power = softpower)
  print('NOTE: Done calculating adjacency values')

  # NOTE
  # To minimize effects of noise and spurious associations we transform adjacency matrix
  # into a Topological Overlap Matrix and calculate dissimilarity.
  # Turn adjacency into topological overlap
  # Fails unless you have an absurd amount of RAM, necessary to move to HPCC at this point
  print('NOTE: Calculating dissTOM now')
  dissTOM <<- calc_diss_TOM(adjacency)
  rm(adjacency)  # clear space
  
  geneTree <<- cluster_using_TOM(dissTOM)
  
  trim_and_color(geneTree, dissTOM)
  
  merging_similar_expression(working_data, dynamicColors, geneTree)
  
  names_and_colors(working_data)
  
  geneInfo_dframe <<- data.frame(Gene_Names = names(working_data),
                                 moduleColor = moduleColors)
  write.table(
    geneInfo_dframe,
    'Genes_and_ModuleColors.tsv',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
  )
   
  writing_intermediate_data = t(working_data)
  writing_intermediate_data = tibble::rownames_to_column(working_data, 'Gene_Names')
  write.table(
    writing_intermediate_data,
    'Intermediate_Data_TPM.tsv',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
  )
  
}

# Choose a soft power of 3 because that was best for the all test set
process(load_all_hap_hpcc(), 'All_Genes', 3)

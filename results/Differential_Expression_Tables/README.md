# How to Interpret Results:
Each sheet (and cell) in the results file is derived from the `_Summary.txt` files that were outputted from EdgeR.
The data in the file is derived from the all haplotype dataset using the FDR correction (BH, not Bonferroni).
For each cell where there is a data point, there is a 3-tuple containing the number of genes that are (down-regulated, not significantly different in expression, and up-regulated).
The genome that is column-major is the genome in which the genes are up or down-regulated. For example, if the genome along the columns is **Liberty**, and there is a cell that says (1560, 91420, 603), then there are 1560 genes down-regulated in the Liberty genome, and 603 up-regulated in the Liberty genome.

# By Scott
# NOTE we did not take GO hierarchy into account because since this dataset is already downstream
# of a lot of filtering, we did not want to penalize the model further. 

# NOTE each installation command must be done separately
# Run these commands if it is the first time
# install.packages("BiocManager")
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")
# install.packages("tidyverse")

renv::restore()
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(Rgraphviz))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tools))

############################################
# OBJECTIVE:
	# Perform a GO enrichment on the DEG data.

# Steps:
	# Gather the Draper DEG data (days 1-7) and join the genes into 1 file TODO check
	# Gather the Liberty DEG data (days 1-7) and join the genes into 1 file TODO check
	# Define the 'my_interesting_genes' variable, this will change according to
		# Draper or Liberty
	# Define the "GO Universe" (master_genes) that you are searching against
	# Perform a single GO enrichment

############################################
# Define input args and do simple set-up

# Check to see you have the args
args= commandArgs(trailingOnly=TRUE)

# Draper-Arabidopsis DEGs
draper = read.table(file=args[1])
liberty = read.table(file=args[2])



# Define your mapping (gene universe), MAGIC input arg 2
geneID2GO = readMappings(file=args[3], sep='\t') # NB global
master_genes = names(geneID2GO)  # NB global

# output directory
output_dir = args[4]
doc_dir = args[5]

############################################
# Begin the formal analysis

function_run_topgo = function(master_genes, geneID2GO, my_interesting_genes, module_name, ontology_group, output_dir){
	geneList = factor(as.integer(master_genes %in% my_interesting_genes))
	names(geneList) = master_genes # NB give it names

	GOdata = new('topGOdata',
		     ontology = ontology_group,
		     allGenes = geneList,
		     annot = annFUN.gene2GO,
		     gene2GO = geneID2GO)

	# NOTE get the number of genes in this interesting list
	# NOTE how many genes am I losing when I perform this step?
	# Do some genes not have the 'MF' ontology or something?
	# print(GOdata)  # NB shows the number of genes lost

	sig_genes = sigGenes(GOdata)  # the significant genes
	num_sig_genes = numSigGenes(GOdata) # NB the no. of signifcant genes,
	# not necessarily the nodes in the graph
	num_nodes = length(usedGO(object=GOdata))  # NB all nodes in the graph,
	# not necessarily the significant nodes
	# FUTURE would this be useful information to store?

	#resultFisher_weight = runTest(GOdata, algorithm='weight01', statistic='fisher')
	resultFisher_classic = runTest(GOdata, algorithm='classic', statistic='fisher')


	# NB from Nolan:
		# algorithm="classic" WILL NOT take GO hierarchy into account
	    	# ^ The limitation of using "classic" is that all genes
       		# annotated to a GO terms will be automatically annotated
       		# to its parents as well, therefore a GO term might look
       		# enriched just because its children are enriched
	    	# algorithm="weight01" WILL take GO hierarchy into account
		# These p-values have not been corrected for multiple testing

	# NB from Pat:
       		# After discussion with Pat, since I already have a cutoff applied
       		# to my modules, and the GO hieracrhy is useful information, I don't want
       		# to penalize the stats further, so I am going with classic and won't do
       		# an FDR.

	# list the top significant GO terms
    	allRes = GenTable(GOdata, classicFisher=resultFisher_classic, 
			  orderBy='classicFisher',
			  ranksOf='classicFisher',
			  topNodes=num_nodes,
			  numChar=1000) # NB avoid truncation of string
			      
			      
	# MAGIC p-val threshold
       	pval_threshold = 0.05
	# NB apply the p-val threshold
	allRes_significant = allRes[which(allRes$classicFisher < pval_threshold),]


	# NB add a string to a new column that says 'yes' if we have
       	# more significant than expected
	allRes_significant$Overrepresented = ifelse(allRes_significant$Significant > allRes_significant$Expected, 'Y', 'N')

	# NB A GO term is synonymoous with node
	# NB node size defaults to 1, no pruning is performed
	# This is the minimum number of genes annotated to a go

	# MAGIC filename
	write.table(allRes_significant, file=file.path(output_dir, paste(module_name, '_', ontology_group, "_",pval_threshold,"_pval_sig_GO_Summary_Table.tsv", sep='')), sep='\t', quote=FALSE, row.names=FALSE)

}	

# Running the functions (needs to be done once for Draper and once for Liberty
function_run_topgo(master_genes, geneID2GO, as.character(draper$V1), 'Draper', 'BP', output_dir)
function_run_topgo(master_genes, geneID2GO, as.character(liberty$V1), 'Liberty', 'BP', output_dir)

# NB this outputs the session information for easy package management.
# Commented out because it is redundant with the other TopGO script
# Perhaps better than outputting a conda env
# MAGIC get nice doc dir
#sink(file=paste(doc_dir, "TopGO_sessionInfo.txt", sep='/'))
#sessionInfo()
#sink()

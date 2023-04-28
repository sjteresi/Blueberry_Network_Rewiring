# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/data
DEV_RESULTS := $(ROOT_DIR)/results
DEV_DOCUMENTATION := $(ROOT_DIR)/doc

# NOTE, this is from the preceding project
DEV_DIFFEXDIR := $(ROOT_DIR)/../Blueberry_RNA_Seq_Expression_Analysis/results/Diff_Ex/EdgeR_Output

DEV_GENE_ANNOTATION := $(DEV_DATA)/V_corymbosum_v1.0_geneModels.gff

# Orthology analysis related paths
DEV_ORTHOLOGY_DATA := $(DEV_DATA)/orthology_data
DEV_SYNTELOGS := $(DEV_ORTHOLOGY_DATA)/SynMap_Arabidopsis_Blueberry_8_12_2020.txt
DEV_HOMOLOGS := $(DEV_ORTHOLOGY_DATA)/At-Blueberry.blast
DEV_SYNTENY_HOMOLOGY_TABLE := $(DEV_ORTHOLOGY_DATA)/Synteny_Homology_Table.tsv 


# FPKM/TPM related paths
DEV_EXPRESSION_OUT_DIR := $(DEV_DATA)/FPKM_TPM
DEV_COLLATED_COUNT_FILE := $(DEV_DATA)/AllCounts_Blueberry.tsv


# WGCNA analysis related paths
DEV_WGCNA_OUT_DIR := $(DEV_DATA)/WGCNA
# Output for later
DEV_WGCNA_GENES_AND_MODULES := $(DEV_WGCNA_OUT_DIR)/Genes_and_ModuleColors.tsv

# Module filtering related paths
DEV_MODULES_IN_AT := $(DEV_DATA)/Modules/modulecolors_AT

# TopGO related paths
DEV_DOWNLOADED_GO_UNIVERSE := $(DEV_DATA)/ATH_GO_GOSLIM.txt
DEV_FILTERED_GO_OUTPUT := $(DEV_DATA)/GO/ArabidopsisGene_w_GO.tsv
DEV_FILTERED_GO_OUTPUT_DIR := $(DEV_DATA)/GO/

DEV_DRAPER_DEG := $(DEV_RESULTS)/DEGs/Unique_and_Shared_DEGs/Unique_Draper_All.txt
DEV_LIBERTY_DEG := $(DEV_RESULTS)/DEGs/Unique_and_Shared_DEGs/Unique_Liberty_All.txt

setup:
	mkdir -p data \
	       	data/Arab_Blueberry_BLAST_DB \
	       	data/orthology_data

# NB, cluster script
.PHONY: blastall
blastall:
	sbatch $(ROOT_DIR)/src/orthology_analysis/blastall.sb	

# NB creates Filtered_Syntelogs.tsv, Filtered_Homologs.tsv, and Synteny_Homology_Table.tsv
.PHONY: gen_ortholog_table
gen_ortholog_table:
	python $(ROOT_DIR)/src/orthology_analysis/filter_orthologs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_ORTHOLOGY_DATA)

.PHONY: gen_FPKM_table
gen_FPKM_table:
	mkdir -p $(DEV_EXPRESSION_OUT_DIR)
	python $(ROOT_DIR)/src/FPKM_TPM/process_fpkm.py $(DEV_GENE_ANNOTATION) $(DEV_COLLATED_COUNT_FILE) $(DEV_EXPRESSION_OUT_DIR)

.PHONY: gen_TPM_table
gen_TPM_table:
	mkdir -p $(DEV_EXPRESSION_OUT_DIR)
	python $(ROOT_DIR)/src/FPKM_TPM/process_tpm.py $(DEV_GENE_ANNOTATION) $(DEV_COLLATED_COUNT_FILE) $(DEV_EXPRESSION_OUT_DIR)

# NOTE WGCNA must be run on cluster. Just putting the command here in the Makefile for reference of order of code run.
# A feature of it being run on the cluster and in R is that some of the code-paths are hard-coded in the script
.PHONY: run_WGCNA
run_WGCNA:
	mkdir -p $(DEV_WGCNA_OUT_DIR)/WGCNA
	sbatch $(ROOT_DIR)/src/WGCNA/run_WGCNA.sb

# Convert modules (blueberry genes) to Arabidopsis genes
# NOTE this step is necessary for TopGO even though the summary table produces a more legible table
.PHONY: blueberry_module_conversion_to_arabidopsis
blueberry_module_conversion_to_arabidopsis:
	mkdir -p $(DEV_DATA)/Modules
	# Note sub dirs are made in the Python script
	python $(ROOT_DIR)/src/modules/filter_modules.py $(DEV_WGCNA_GENES_AND_MODULES) $(DEV_SYNTENY_HOMOLOGY_TABLE) $(DEV_DATA)/Modules

# Work on GO:
# Distill the raw gene universe GO file down into a format for TopGO
.PHONY: filter_GO
filter_GO:
	mkdir -p $(DEV_FILTERED_GO_OUTPUT_DIR)
	python $(ROOT_DIR)/src/TopGO/generate_gene_w_GO_term.py $(DEV_DOWNLOADED_GO_UNIVERSE) $(DEV_FILTERED_GO_OUTPUT_DIR)

# R code on personal
.PHONY: topGO_modules
topGO_modules:
	mkdir -p $(DEV_DATA)/GO/TopGO_Modules
	Rscript $(ROOT_DIR)/src/TopGO/topGO_blueberry.R $(DEV_MODULES_IN_AT) $(DEV_FILTERED_GO_OUTPUT) $(DEV_DATA)/GO/TopGO_Modules $(DEV_DOCUMENTATION)

.PHONY: filter_proteins
filter_proteins:
	mkdir -p $(DEV_DATA)/proteins
	python $(ROOT_DIR)/src/proteins/protein_table.py $(DEV_DATA)/proteins/TAIR10_pep_20101214.txt $(DEV_DATA)/proteins

# Gene Stats
.PHONY: gene_stats_all
gene_stats_all:
	python $(ROOT_DIR)/src/orthology_analysis/identification_proportions.py $(DEV_GENE_ANNOTATION) $(DEV_SYNTENY_HOMOLOGY_TABLE)

.PHONY: generate_exp_table_melanie
generate_exp_table_melanie:
	python3 $(ROOT_DIR)/src/exp_table_melanie.py $(DEV_DATA)/FPKM_TPM/Blueberry_FPKM_All_Haplotype.tsv $(DEV_ORTHOLOGY_DATA)/Synteny_Homology_Table.tsv $(DEV_ORTHOLOGY_DATA)


# NB this is being pipelined into another project
# Filter Log2FC by modules
.PHONY: module_log2fc_overlap
module_log2fc_overlap:
	mkdir -p $(DEV_DATA)/module_overlap/module_log2fc_overlap/
	python $(ROOT_DIR)/src/module_overlap/module_log2fc_overlap.py \
	       $(DEV_DATA)/Log_2FC_Melanie/ \
	       $(DEV_WGCNA_GENES_AND_MODULES) \
	       $(DEV_DATA)/module_overlap/module_log2fc_overlap/


# Filter DEGs by Modules, the DEG data was generated in the preceding project
.PHONY: module_deg_overlap
module_deg_overlap:
	mkdir -p $(DEV_DATA)/module_overlap/module_deg_overlap/
	python $(ROOT_DIR)/src/module_overlap/module_deg_overlap.py \
	       $(DEV_DIFFEXDIR)/All_Hap/FDR/ \
	       $(DEV_WGCNA_GENES_AND_MODULES) \
	       $(DEV_SYNTENY_HOMOLOGY_TABLE) \
	       $(DEV_DATA)/module_overlap/module_deg_overlap/

# Make figures and tables relating to DEG representation at various time points
# NOTE, these graphs were made and Melanie "re-did" them in illustrator to touch them up
.PHONY: deg_time_points
deg_time_points:
	mkdir -p $(DEV_DATA)/DEGs/Unique_and_Shared_DEGs
	python $(ROOT_DIR)/src/DEG_Analysis/deg_time_points.py \
	       $(DEV_DIFFEXDIR)/All_Hap/FDR/ \
	       $(DEV_DATA)/GO/GO_ID_w_Term.tsv \
	       $(DEV_SYNTENY_HOMOLOGY_TABLE) \
	       $(DEV_DATA)/DEGs/

.PHONY: module_expression_graphs
module_expression_graphs:
	mkdir -p $(DEV_DATA)/module_expression/
	python $(ROOT_DIR)/src/modules/module_expression_graphs.py \
	$(DEV_WGCNA_GENES_AND_MODULES) \
	$(DEV_DATA)/FPKM_TPM/Blueberry_TPM_All_Haplotype.tsv \
	$(DEV_DATA)/module_expression/

.PHONY: deg_qtl
deg_qtl:
	mkdir -p $(DEV_DATA)/QTL
	python $(ROOT_DIR)/src/QTL/deg_qtl.py \
	       $(DEV_DIFFEXDIR)/All_Hap/FDR/ \
	       $(DEV_DATA)/QTL/QTL_genes_of_interest.csv \
	       $(DEV_DATA)/module_expression/Mean_Expression_TPM.tsv \
	       $(DEV_SYNTENY_HOMOLOGY_TABLE) \
	       $(DEV_DATA)/Log_2FC_Melanie/ \
	       $(DEV_DATA)/proteins/Filtered_Arabidopsis_Protein_Info.tsv \
	       $(DEV_DATA)/GO/GO_ID_w_Term.tsv \
	       $(DEV_DATA)/QTL/

# NOTE I hate R and its package management is horrible. You cannot execute this 
# Makefile command. You must be in the directory of the R script and execute it
# yourself. If you don't do it that way, it won't load the R environment of installed
# packages correctly. SO, that means this Makefile command is mainly just here for
# reference.
# I am also runnning this on the HPCC, so to load R itself:
# module load GCC/8.3.0 OpenMPI/3.1.4 R/3.6.2
topGO_deg:
	mkdir -p $(DEV_DATA)/GO/TopGO_DEGs
	Rscript $(ROOT_DIR)/src/TopGO/topGO_DEGs.R $(DEV_DRAPER_DEG) $(DEV_LIBERTY_DEG) $(DEV_FILTERED_GO_OUTPUT) $(DEV_DATA)/GO/TopGO_DEGs $(DEV_DOCUMENTATION)

# Filter GO by Modules
# NB the List_Top* files were given by Melanie
.PHONY: module_go_overlap
module_go_overlap:
	mkdir -p $(DEV_DATA)/module_overlap/module_go_overlap/
	python $(ROOT_DIR)/src/module_overlap/module_go_overlap.py \
	       $(DEV_DATA)/GO/TopGO_Modules \
	       $(DEV_DATA)/module_overlap/module_go_overlap/List_Top_10_GO.tsv \
	       $(DEV_DATA)/module_overlap/module_go_overlap/List_Top_11_Modules.tsv \
	       $(DEV_DATA)/Log_2FC_Melanie/ \
	       $(DEV_SYNTENY_HOMOLOGY_TABLE) \
	       $(DEV_DATA)/GO/GO_ID_w_Term.tsv \
	       $(DEV_DATA)/proteins/Filtered_Arabidopsis_Protein_Info.tsv \
	       $(DEV_WGCNA_GENES_AND_MODULES) \
	       $(DEV_DATA)/module_overlap/module_go_overlap/


# NB NOT IMPLEMENTED
# CANDIDATE FOR DELETION
create_module_table:
	$(ROOT_DIR)/src/modules/make_module_deg_table.py $(DEV_WGCNA_GENES_AND_MODULES) $(DEV_DIFFEXDIR)/ "All_Hap" "FDR" $(DEV_RESULTS)/Modules

# NB NOT IMPLEMENTED
# CANDIDATE FOR DELETION
create_Log_2FC_histogram:
	$(ROOT_DIR)/src/modules/log2fc_hist.py $(DEV_RESULTS)/Log_2FC_Melanie/Melanie_Log_2FC_Filtered.tsv $(DEV_RESULTS)/Log_2FC_Melanie/histograms.png


get_lists_of_Arabidopsis_degs:
	awk '{if($$1 != "No_Ortholog" && $$1 != "Arabidopsis_Gene") print $$1}' $(ROOT_DIR)/results/DEGs/Unique_and_Shared_DEGs/Unique_Draper_* > $(ROOT_DIR)/results/DEGs/Unique_and_Shared_DEGs/Unique_Draper_All.txt
	awk '{if($$1 != "No_Ortholog" && $$1 != "Arabidopsis_Gene") print $$1}' $(ROOT_DIR)/results/DEGs/Unique_and_Shared_DEGs/Unique_Liberty* > $(ROOT_DIR)/results/DEGs/Unique_and_Shared_DEGs/Unique_Liberty_All.txt





sync_local_to_remote_data:
	rsync -ave ssh /home/scott/Documents/Uni/Research/Projects/Blueberry_Network_Rewiring/data --chmod=Dg+s teresisc@rsync.hpcc.msu.edu:/mnt/research/edgerpat_lab/Scotty/Blueberry_Network_Rewiring/

sync_local_to_remote_results:
	rsync -ave ssh /home/scott/Documents/Uni/Research/Projects/Blueberry_Network_Rewiring/results --chmod=Dg+s teresisc@rsync.hpcc.msu.edu:/mnt/research/edgerpat_lab/Scotty/Blueberry_Network_Rewiring/

sync_hpcc_to_onedrive:
	# MUST be standing in root folder for project
	ml Rclone
	rclone sync . remote:HPCC_Mirror/Blueberry_Network_Rewiring/ --exclude=renv/** --exclude=.git/** -P


# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/data
DEV_RESULTS := $(ROOT_DIR)/results
DEV_DIFFEXDIR := $(DEV_DATA)/Diff_Ex/EdgeR_Output

DEV_GENE_ANNOTATION := $(DEV_DATA)/V_corymbosum_v1.0_geneModels.gff

# Orthology analysis related paths
DEV_SYNTELOGS := $(DEV_DATA)/synmap_out_8_12_2020.txt
DEV_HOMOLOGS := $(DEV_DATA)/At-Blueberry.blast
DEV_ORTHOLOGY_OUT_DIR := $(DEV_RESULTS)/Arabidopsis_Blueberry_Orthology
# Output for later
DEV_SYNTENY_HOMOLOGY_TABLE := $(DEV_ORTHOLOGY_OUT_DIR)/Synteny_Homology_Table.tsv 


# FPKM/TPM related paths
DEV_EXPRESSION_OUT_DIR := $(DEV_RESULTS)/FPKM_TPM
DEV_COLLATED_COUNT_FILE := $(DEV_DATA)/AllCounts_Blueberry.tsv


# WGCNA analysis related paths
DEV_WGCNA_OUT_DIR := $(DEV_RESULTS)/WGCNA
# Output for later
DEV_WGCNA_GENES_AND_MODULES := $(DEV_WGCNA_OUT_DIR)/Genes_and_ModuleColors.tsv

# TopGO related paths
DEV_DOWNLOADED_GO_UNIVERSE := $(DEV_DATA)/ATH_GO_GOSLIM.txt
DEV_FILTERED_GO_OUTPUT := $(DEV_RESULTS)/GO/ArabidopsisGene_w_GO.tsv


.PHONY: dev help

gen_ortholog_table:
	mkdir -p $(DEV_ORTHOLOGY_OUT_DIR)
	$(ROOT_DIR)/src/Arabidopsis_Blueberry_Orthology/filter_orthologs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_ORTHOLOGY_OUT_DIR)

gen_FPKM_table:
	mkdir -p $(DEV_EXPRESSION_OUT_DIR)
	python $(ROOT_DIR)/src/FPKM_TPM/process_fpkm.py $(DEV_GENE_ANNOTATION) $(DEV_COLLATED_COUNT_FILE) $(DEV_EXPRESSION_OUT_DIR)

gen_TPM_table:
	mkdir -p $(DEV_EXPRESSION_OUT_DIR)
	python $(ROOT_DIR)/src/FPKM_TPM/process_tpm.py $(DEV_GENE_ANNOTATION) $(DEV_COLLATED_COUNT_FILE) $(DEV_EXPRESSION_OUT_DIR)

# NOTE WGCNA must be run on cluster. Just putting the command here in the Makefile for reference of order of code run.
# A feature of it being run on the cluster and in R is that some of the code-paths are hard-coded in the script
run_WGCNA:
	mkdir -p $(DEV_RESULTS)/WGCNA
	sbatch $(ROOT_DIR)/src/WGCNA/run_WGCNA.sb

# Work on GO:
# Distill the raw gene universe GO file down into a format for TopGO
GO:
	mkdir -p $(DEV_RESULTS)/GO
	python $(ROOT_DIR)/src/TopGO/generate_gene_w_GO_term.py $(DEV_DOWNLOADED_GO_UNIVERSE) $(DEV_FILTERED_GO_OUTPUT)

# Run TopGO to get over/under representation of terms
topGO:
	Rscript $(ROOT_DIR)/src/TopGO/topGO_blueberry.R $(DEV_DATA)/WGCNA_Data/modulecolors_AT $(DEV_RESULTS)/GO/ArabidopsisGene_w_GO.tsv $(DEV_RESULTS)/GO

#----------------------------------------#
# Master summary table
summary_table:
	mkdir -p $(DEV_RESULTS)/Summary_Diff_Ex_Modules
	python $(ROOT_DIR)/src/summary_table.py \
	       $(DEV_DATA)/WGCNA_Data/Final_WGCNA/Genes_and_ModuleColors.tsv \
	       $(DEV_DATA)/AtBB/data_output/Synteny_Homology_Table.tsv \
	       $(DEV_DATA)/Diff_Ex/EdgeR_Output/All_Hap/FDR/ \
	       $(DEV_RESULTS)/GO/ArabidopsisGene_w_GO.tsv \
	       $(DEV_DATA)/TPM/Blueberry_TPM_All.tsv \
	       $(DEV_RESULTS)/Complete_Gene_Summary_Table.tsv \
	       $(DEV_RESULTS)/Summary_Diff_Ex_Modules


#----------------------------------------#
# Gene Stats
gene_stats_all:    # execute gene stats
	$(ROOT_DIR)/src/gene_stats/operations.py $(DEV_GENOME) $(DEV_ORTHOLOGY)

#----------------------------------------#

generate_exp_table_melanie:
	mkdir -p $(DEV_RESULTS)/melanie_exp_syn_table
	python3 $(ROOT_DIR)/src/exp_table_melanie.py $(DEV_DATA)/FPKM/Blueberry_FPKM_All.tsv $(DEV_DATA)/AtBB/data_output/merged_homo_and_syn.tsv $(DEV_RESULTS)/melanie_exp_syn_table



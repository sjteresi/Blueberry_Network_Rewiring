# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/data
DEV_RESULTS := $(ROOT_DIR)/results
DEV_DIFFEXDIR := $(DEV_DATA)/Diff_Ex/EdgeR_Output

# Orthology analysis related paths
DEV_SYNTELOGS := $(DEV_DATA)/synmap_out_8_12_2020.txt
DEV_HOMOLOGS := $(DEV_DATA)/At-Blueberry.blast
DEV_ORTHOLOGY_OUT_DIR := $(DEV_RESULTS)/Arabidopsis_Blueberry_Orthology

DEV_GENOME := $(DEV_DATA)/Genome/V_corymbosum_v1.0_geneModels.gff
DEV_ORTHOLOGY := $(DEV_DATA)/AtBB/data_output/merged_homo_and_syn.tsv

# DEV_RESULTS := $(DEV_DATA)/results

.PHONY: dev help

gen_ortholog_table:
	mkdir -p $(DEV_ORTHOLOGY_OUT_DIR)
	$(ROOT_DIR)/src/Arabidopsis_Blueberry_Orthology/filter_orthologs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_ORTHOLOGY_OUT_DIR)

# Convert modules (blueberry genes) to Arabidopsis genes
module_conversion:
	mkdir -p $(DEV_RESULTS)/Modules
	python $(ROOT_DIR)/src/modules/filter_modules.py $(DEV_DATA)/WGCNA_Data/Final_WGCNA/Genes_and_ModuleColors.tsv $(DEV_DATA)/AtBB/data_output/Synteny_Homology_Table.tsv $(DEV_RESULTS)/Modules

create_diff_ex_tables:
	mkdir -p $(DEV_RESULTS)/Differential_Expression_Tables/
	python $(ROOT_DIR)/src/diff_ex_tables/diff_ex_tables.py $(DEV_DATA)/Diff_Ex/EdgeR_Output/All_Hap/FDR/ $(DEV_RESULTS)/Differential_Expression_Tables


# Work on GO:
# Distill the raw gene universe /  GO file down into a format for TopGO
GO:
	mkdir -p $(DEV_RESULTS)/GO
	python $(ROOT_DIR)/src/TopGO/generate_gene_w_GO_term.py $(DEV_DATA)/GO/ATH_GO_GOSLIM.txt $(DEV_RESULTS)/GO

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

# TODO add methods for FPKM


fpkm:                ## execute fpkm code with our data


generate_exp_table_melanie:
	mkdir -p $(DEV_RESULTS)/melanie_exp_syn_table
	python3 $(ROOT_DIR)/src/exp_table_melanie.py $(DEV_DATA)/FPKM/Blueberry_FPKM_All.tsv $(DEV_DATA)/AtBB/data_output/merged_homo_and_syn.tsv $(DEV_RESULTS)/melanie_exp_syn_table


help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

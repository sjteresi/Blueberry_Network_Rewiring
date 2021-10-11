# scripts for development
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/../Blueberry_Data
DEV_RESULTS := $(ROOT_DIR)/results
DEV_SYNTELOGS := $(DEV_DATA)/AtBB/data_input/synmap_out_8_12_2020.txt
DEV_HOMOLOGS := $(DEV_DATA)/AtBB/data_input/At-Blueberry.blast
DEV_DIFFEXDIR := $(DEV_DATA)/Diff_Ex/EdgeR_Output

DEV_GENOME := $(DEV_DATA)/Genome/V_corymbosum_v1.0_geneModels.gff
DEV_ORTHOLOGY := $(DEV_DATA)/AtBB/data_output/merged_homo_and_syn.tsv

# DEV_RESULTS := $(DEV_DATA)/results

.PHONY: dev help

# All Hap ATBB
at_bb_allbf:                ## execute orthology code with our data
	$(ROOT_DIR)/scripts/At_BB/generate_pairs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_DIFFEXDIR)/All_Hap/Bonferroni Bonferroni All

at_bb_allfdr:                ## execute orthology code with our data
	$(ROOT_DIR)/scripts/At_BB/generate_pairs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_DIFFEXDIR)/All_Hap/FDR FDR All

# Single Hap ATBB
at_bb_singlebf:                ## execute orthology code with our data
	$(ROOT_DIR)/scripts/At_BB/generate_pairs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_DIFFEXDIR)/Single_Hap/Bonferroni Bonferroni Single
at_bb_singlefdr:                ## execute orthology code with our data
	$(ROOT_DIR)/scripts/At_BB/generate_pairs.py $(DEV_SYNTELOGS) $(DEV_HOMOLOGS) $(DEV_DIFFEXDIR)/Single_Hap/FDR FDR Single

# Convert modules (blueberry genes) to Arabidopsis genes
module_conversion:
	mkdir -p $(DEV_RESULTS)/Modules
	python $(ROOT_DIR)/scripts/modules/filter_modules.py $(DEV_DATA)/WGCNA_Data/Final_WGCNA/Genes_and_ModuleColors.tsv $(DEV_DATA)/AtBB/data_output/Synteny_Homology_Table.tsv $(DEV_RESULTS)/Modules

# Gene Stats
gene_stats_all:    # execute gene stats
	$(ROOT_DIR)/scripts/gene_stats/operations.py $(DEV_GENOME) $(DEV_ORTHOLOGY)

#----------------------------------------#

# TODO add methods for FPKM


fpkm:                ## execute fpkm code with our data


generate_exp_table_melanie:
	mkdir -p $(DEV_RESULTS)/melanie_exp_syn_table
	python3 $(ROOT_DIR)/scripts/exp_table_melanie.py $(DEV_DATA)/FPKM/Blueberry_FPKM_All.tsv $(DEV_DATA)/AtBB/data_output/merged_homo_and_syn.tsv $(DEV_RESULTS)/melanie_exp_syn_table


help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

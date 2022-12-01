#!/usr/bin/env python3

"""
Takes the list of interesting genes from the results of the QTL analysis and
determines which genes are also DEGs
"""

__author__ = "Scott Teresi"

# General imports
import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

# Import reader functions
from src.read_tables_and_dir import (
    read_DEG_dir,
    read_all_DEG_direction_files,
    read_GO_ID_w_term,
    read_FPKM_or_TPM,
    read_synteny_homology_table,
    read_log_2fc_dir,
    read_all_log_2fc_files,
    read_protein_table,
    convert_direction_integer_to_string,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    # NB input args correspond to Makefile, must be ordered correctly
    parser.add_argument("DEG_dir", type=str, help="output from EdgeR")
    parser.add_argument("qtl_genes", type=str, help="output from qtl anaysis")
    parser.add_argument("expression_matrix", type=str, help="expression file")
    parser.add_argument("orthologs", type=str, help="synteny homology table")
    parser.add_argument("log2fc_dir", type=str, help="log_2fc_info")
    parser.add_argument("protein_info", type=str, help="protein_info")
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.DEG_dir = os.path.abspath(args.DEG_dir)
    args.qtl_genes = os.path.abspath(args.qtl_genes)
    args.expression_matrix = os.path.abspath(args.expression_matrix)
    args.orthologs = os.path.abspath(args.orthologs)
    args.log2fc_dir = os.path.abspath(args.log2fc_dir)
    args.protein_info = os.path.abspath(args.protein_info)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE begin reading data
    list_of_deg_paths = read_DEG_dir(args.DEG_dir)
    list_of_deg_pandas, analysis_contexts = read_all_DEG_direction_files(
        list_of_deg_paths
    )
    list_of_deg_pandas_vals_as_str = list(
        map(convert_direction_integer_to_string, list_of_deg_pandas, analysis_contexts)
    )

    # Table of blueberry genes and the direction of differential regulation per
    # context
    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )

    # Table of blueberry genes and their Arabidopsis orthologs
    orthologs = read_synteny_homology_table(args.orthologs)

    # Table of Log2FC info
    log2fc_panda = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        read_all_log_2fc_files(read_log_2fc_dir(args.log2fc_dir)),
    )
    log2fc_col = [col for col in log2fc_panda.columns if "FC" in col]

    # Get protein information
    protein_table = read_protein_table(args.protein_info)

    # QTL genes
    qtl_genes = pd.read_csv(args.qtl_genes, header="infer", sep=",")
    print(
        f"This is how many unique blueberry QTL genes were identified by collaborators: {qtl_genes['Blueberry_Gene'].nunique()}"
    )

    # Merge in Log2FC info, the ortholog info, and the protein info
    qtl_genes = qtl_genes.merge(log2fc_panda, on="Blueberry_Gene", how="left")
    qtl_genes = qtl_genes.merge(orthologs, on="Blueberry_Gene", how="left")
    qtl_genes = qtl_genes.merge(protein_table, on="Arabidopsis_Gene", how="left")
    qtl_genes.fillna("NA", inplace=True)

    # Get number of unique blueberry genes that have an AT ortholog
    booled = qtl_genes["Arabidopsis_Gene"] != "NA"
    subsetted = qtl_genes[booled]
    num_unique_bb_genes_with_at_ortho = subsetted["Blueberry_Gene"].nunique()
    print(
        f"This is how many unique blueberry genes have an AT ortholog: {num_unique_bb_genes_with_at_ortho}"
    )

    # Get number of unique blueberry genes that are found in at least one
    # log2fc context
    booled = (qtl_genes[log2fc_col] != "NA").any(axis=1)
    subsetted = qtl_genes[booled]
    num_unique_bb_genes_with_at_least_log2fc = subsetted["Blueberry_Gene"].nunique()
    print(
        f"This is how many unique blueberry genes were log2fc expressed in at least one context: {num_unique_bb_genes_with_at_least_log2fc}"
    )

    # Get general expression info
    expression_data = read_FPKM_or_TPM(args.expression_matrix)
    expression_data.set_index("Blueberry_Gene", inplace=True)
    columns_of_interest = [
        col for col in expression_data.columns.to_list() if "All_" not in col
    ]
    any_expression_bool = (expression_data[columns_of_interest] != 0).any(axis=1)
    any_expression_at_all = expression_data[any_expression_bool]
    qtl_genes = pd.merge(
        any_expression_at_all, qtl_genes, on="Blueberry_Gene", how="inner"
    )
    qtl_genes_with_any_expression_at_all = qtl_genes.copy(deep=True)
    qtl_genes_with_any_expression_at_all.drop(columns=log2fc_col, inplace=True)
    qtl_genes_with_any_expression_at_all.sort_index(axis=1, inplace=True)
    num_unique_bb_genes_with_any_expression = qtl_genes_with_any_expression_at_all[
        "Blueberry_Gene"
    ].nunique()

    print(
        f"This is how many unique blueberry genes were generally expressed in at least one context: {num_unique_bb_genes_with_any_expression}"
    )

    print(qtl_genes)
    # Save the data
    qtl_genes.sort_index(axis=1, inplace=True)
    qtl_genes.to_csv(
        os.path.join(args.output_dir, "candidate_genes_qtl_table.tsv"),
        sep="\t",
        header=True,
        index=False,
    )

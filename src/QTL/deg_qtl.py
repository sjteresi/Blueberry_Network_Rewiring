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
    read_synteny_homology_table,
    convert_direction_integer_to_string,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    # NB input args correspond to Makefile, must be ordered correctly
    parser.add_argument("DEG_dir", type=str, help="output from EdgeR")
    parser.add_argument("qtl_genes", type=str, help="output from qtl anaysis")
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.DEG_dir = os.path.abspath(args.DEG_dir)
    args.qtl_genes = os.path.abspath(args.qtl_genes)
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
    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )

    qtl_genes = pd.read_csv(args.qtl_genes, header="infer", sep=",")

    qtl_genes = pd.merge(
        merged_deg_no_module_ID, qtl_genes, on="Blueberry_Gene", how="inner"
    )

    qtl_genes.set_index("Blueberry_Gene", inplace=True)
    qtl_genes.to_csv(
        os.path.join(args.output_dir, "all_qtl_genes.tsv"),
        sep="\t",
        header=True,
        index=True,
    )
    non_scaffold_col = [col for col in qtl_genes.columns if "Scaffold_ID" not in col]
    # Begin filtering pandaframe by whether or not we observe a gene being a
    # DEG in at least once comparison context
    boolean_filtered_degs = qtl_genes[non_scaffold_col] != "No_Change"
    bool_mask = boolean_filtered_degs.any(axis=1)
    at_least_one_deg = qtl_genes[bool_mask]
    at_least_one_deg.to_csv(
        os.path.join(args.output_dir, "at_least_one_deg_qtl_genes.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

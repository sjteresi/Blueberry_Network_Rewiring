#!/usr/bin/env python3

"""
Filter DEG representation data. Determine overlap with modules from
WGCNA.

What modules have the highest or lowest amount of DEGs?
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

from src.module_overlap.module_log2fc_overlap import (
    calc_gene_counts_per_module,
    top_n_modules_and_percent_per_context,
)

from src.read_tables_and_dir import (
    read_DEG_dir,
    read_all_DEG_direction_files,
    read_synteny_homology_table,
    read_gene_modules_table,
)


def count_direction_of_regulation_per_module_deg(all_merged):
    """
    Converts the Log2FC columns of all_merged into counts per module. Index is
    now the module ID not the gene. There is a count of genes in each context

    Args:
        all_merged (pandas.DataFrame):
            Index:
                RangeIndex:
            Columns:
                Blueberry_Gene: Object, dtype: object
                Number_of_Blueberry_Genes_in_Module: Integer, dtype: int64
                Module_Color: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

    Returns:
        counts (pandas.DataFrame):
            Index:
                Module_Color: Length 89: Object, dtype: object

            Columns:
                Expression Contexts (MAGIC VARIES): Object: dtype: object
    """
    unwanted_columns = [
        "Blueberry_Gene",
        "Module_Color",
        "Number_of_Blueberry_Genes_in_Module",
    ]
    columns_of_direction = all_merged.columns.to_list()
    columns_of_direction = [
        column_name
        for column_name in columns_of_direction
        if column_name not in unwanted_columns
    ]

    counts_per_module_comparison = []
    for direction_column in columns_of_direction:
        x = (
            all_merged.groupby(by="Module_Color", dropna=True)[direction_column]
            .value_counts()
            .unstack()
        )
        comparison_name = x.columns.name
        x = x.add_prefix(comparison_name + "_")
        x = x.rename_axis(None, axis=1)
        counts_per_module_comparison.append(x)
        col_list = [col for col in x.columns if "No_Change" not in col]
        x[direction_column + "_Total_DEGS"] = x[col_list].sum(axis=1)

    # Perform a pandas merge on a list
    counts = reduce(
        lambda left, right: pd.merge(left, right, on="Module_Color", how="outer"),
        counts_per_module_comparison,
    )
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("DEG_dir", type=str, help="output from EdgeR")
    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )

    parser.add_argument(
        "ortholog_pairs",
        type=str,
        help="path to synteny/homology file containing AT BB pairs",
    )
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.DEG_dir = os.path.abspath(args.DEG_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.ortholog_pairs = os.path.abspath(args.ortholog_pairs)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE read the input datasets
    blueberry_genes_and_module_colors_table = read_gene_modules_table(
        args.genes_and_module_colors
    )
    gene_counts_per_module = calc_gene_counts_per_module(
        blueberry_genes_and_module_colors_table
    )
    synteny_homology_table = read_synteny_homology_table(args.ortholog_pairs)

    # NOTE begin the analysis
    list_of_deg_paths = read_DEG_dir(args.DEG_dir)
    list_of_deg_pandas = read_all_DEG_direction_files(list_of_deg_paths)

    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )

    merged_deg_w_module_ID = pd.merge(
        merged_deg_no_module_ID,
        blueberry_genes_and_module_colors_table,
        on="Blueberry_Gene",
        how="outer",
    )
    degs_and_modules = merged_deg_w_module_ID.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )

    modules_and_count_per_context = count_direction_of_regulation_per_module_deg(
        degs_and_modules,
    )

    # Fill the cells with NaN as 0 since we didn't get any genes for that given
    # module and expression context
    modules_and_count_per_context.fillna(0, inplace=True)

    # Drop the 'No_Change' columns
    list_of_no_change_cols = [
        f for f in modules_and_count_per_context.columns if "No_Change" in f
    ]
    modules_and_count_per_context.drop(columns=list_of_no_change_cols, inplace=True)

    # Add counts per module to the table
    modules_and_count_per_context = modules_and_count_per_context.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )
    modules_and_count_per_context.set_index("Module_Color", inplace=True)
    modules_and_count_per_context.sort_index(axis=1, inplace=True, ascending=False)
    # Save the table of the raw counts before we do the division
    modules_and_count_per_context.to_csv(
        os.path.join(args.output_dir, "DEG_counts_in_modules.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

    # Divide to get the percentage
    percentages = modules_and_count_per_context.apply(
        lambda row: row / row["Number_of_Blueberry_Genes_in_Module"], axis=1
    )  # NB calculate percent, make it a new pandas dataframe
    percentages.drop(
        columns=["Number_of_Blueberry_Genes_in_Module"], inplace=True
    )  # NB remove column
    percentages = percentages.add_prefix("Percent_")  # NB add prefix to column
    # NOTE FUTURE
    # Get gene counts per module again because the number of blueberry genes
    # per module column was set to 1 because of division, in future figure out
    # how to do if else statements in pandas apply lambda. Not sure how to keep
    # it from acting on a specific column, would I first need to provide a
    # mask instead of doing if else?
    percentages = percentages.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )

    percentages.set_index("Module_Color", inplace=True)
    percentages.sort_index(axis=1, inplace=True, ascending=False)
    top_n_modules_and_percent_per_context(percentages, args.output_dir)

    # MAGIC filename
    percentages.to_csv(
        os.path.join(args.output_dir, "DEG_percentages_in_modules.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

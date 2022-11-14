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
    convert_direction_integer_to_string,
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

        # NB edge-case where we have zero differentially expressed genes,
        # wouldn't work with the groupby, valuecounts, unstack combo above,
        # resulting in situation where you never get that column. For example,
        # LIB_C2_vs_LIB_T2_Direction_Up doesn't have any genes for this
        # category, which would result in a dataframe without that column,
        # which was confusing. I want the column even if it is full of zeroes
        if "Up" not in x.columns:
            x["Up"] = np.nan
        if "Down" not in x.columns:
            x["Down"] = np.nan

        comparison_name = x.columns.name
        x = x.add_prefix(comparison_name + "_")
        x = x.rename_axis(None, axis=1)
        counts_per_module_comparison.append(x)
        col_list = [col for col in x.columns if "No_Change" not in col]
        x[direction_column + "_Total_DEGs"] = x[col_list].sum(axis=1)

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
    list_of_deg_pandas, analysis_contexts = read_all_DEG_direction_files(
        list_of_deg_paths
    )
    # ----------------------
    # NOTE FUTURE, gettig this 'to_sum' and 'summed' object could be placed in
    # a function for enhanced legibility.
    to_sum = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )
    to_sum.set_index("Blueberry_Gene", inplace=True)
    all_col = to_sum.columns
    for col in all_col:
        to_sum.loc[to_sum[col] == -1, col + "_Direction_Down"] = 1
        to_sum.loc[to_sum[col] == 1, col + "_Direction_Up"] = 1

        # NB, adding the statement for 0 values because in some edge cases
        # where there were zero DEGs, it would never make the column.
        # E.g LIB_C2_vs_LIB_T2_Direction_Up would never get made because there
        # were zero upregulated genes
        to_sum.loc[to_sum[col] == 0, col + "_Direction_Down"] = 0
        to_sum.loc[to_sum[col] == 0, col + "_Direction_Up"] = 0
        to_sum.fillna(0, inplace=True)
        to_sum[col + "_Direction_Total_DEGs"] = (
            to_sum[col + "_Direction_Down"] + to_sum[col + "_Direction_Up"]
        )
        to_sum.drop(columns=[col], inplace=True)

    summed = to_sum.sum(axis=0)

    list_of_deg_pandas_vals_as_str = [
        convert_direction_integer_to_string(panda, context)
        for panda, context in zip(list_of_deg_pandas, analysis_contexts)
    ]

    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas_vals_as_str,
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
    summed = summed.add_prefix("Percent_")
    summed_frame = summed.to_frame("Gene_Counts_Per_Context")

    percentages.set_index("Module_Color", inplace=True)
    transposed_percentages = percentages.T
    percentages_and_counts = transposed_percentages.merge(
        summed_frame, left_index=True, right_index=True, how="outer"
    )
    percentages_and_counts.to_csv(
        os.path.join(args.output_dir, "DEG_percentages_in_modules.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

    # NOTE this is old code but still useful, accessory data tables to
    # 'percentages_and_counts'
    percentages.sort_index(axis=1, inplace=True, ascending=False)
    top_n_modules_and_percent_per_context(percentages, args.output_dir)

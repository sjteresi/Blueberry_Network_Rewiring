#!/usr/bin/env python3

"""
Filter Log2FC gene representation data. Determine overlap with modules from
WGCNA.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

from src.read_tables_and_dir import (
    read_gene_modules_table,
    read_log_2fc_dir,
    read_all_log_2fc_files,
)


def count_direction_of_regulation_per_module_log2fc(all_merged):
    """
    Converts the Log2FC columns of all_merged into counts per module. Index is
    now the module ID not the gene.
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
        # TODO set a better variable name
        x = (
            all_merged.groupby(by="Module_Color", dropna=True)[direction_column]
            .value_counts()
            .unstack()
        )
        comparison_name = x.columns.name
        x = x.add_prefix(comparison_name + "_")
        x = x.rename_axis(None, axis=1)
        counts_per_module_comparison.append(x)

    # Perform a pandas merge on a list
    y = reduce(
        lambda left, right: pd.merge(left, right, on="Module_Color", how="outer"),
        counts_per_module_comparison,
    )
    # TODO set a better variable name
    return y


def calc_gene_counts_per_module(blueberry_genes_and_module_colors_table):
    """
    Takes the filtered output of WGCNA (blueberry genes and their module
    groups), which has already been read in as a pandas frame, and counts the
    number of genes in each module.

    Args:
        blueberry_genes_and_module_colors_table (pandas.DataFrame): Pandas
        dataframe of output from WGCNA. See read_gene_modules_table for more
        detail

    Returns:
        table_of_counts (pandas.DataFrame): Column of 'Module_Color' which is
        string, and column of 'Number_of_Blueberry_Genes_in_Module' as integer.
    """
    table_of_counts = blueberry_genes_and_module_colors_table.groupby(
        "Module_Color"
    ).count()
    table_of_counts.rename(
        columns={"Blueberry_Gene": "Number_of_Blueberry_Genes_in_Module"}, inplace=True
    )
    table_of_counts.reset_index(inplace=True)
    return table_of_counts


def top_n_modules_and_percent_per_context(final, output_dir, n=5):
    """
    TODO fill in
    Output top 5 module ID's for percent genes recovered
    """
    percent_columns = [f for f in final.columns if "Number" not in f]
    for percent_col in percent_columns:
        top_five = final.nlargest(n, columns=[percent_col])
        top_five = top_five.loc[:, top_five.columns.isin([percent_col])]
        top_five.to_csv(
            os.path.join(output_dir, str(percent_col + ".tsv")),
            header=True,
            index=True,
            sep="\t",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # NOTE files will be iterated over
    parser.add_argument(
        "log_2fc_dir", type=str, help="received from Melanie on Google Drive"
    )

    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="output_dir",
    )
    # NB set args
    args = parser.parse_args()
    args.log_2fc_dir = os.path.abspath(args.log_2fc_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    output_dir = os.path.abspath(args.output_dir)

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
    list_log_2fc_files = read_log_2fc_dir(args.log_2fc_dir)
    pandas_log_2fc_files = read_all_log_2fc_files(list_log_2fc_files)

    # NOTE begin the analysis
    # NB do a pandas merge over a list, merge all the log2fc files (which have
    # different columns), but a common blueberry gene column
    log_merged_float_format = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        pandas_log_2fc_files,
    )

    # NB list of expression contexts (important columns) for help later
    expression_contexts = log_merged_float_format.columns.to_list()
    expression_contexts.remove("Blueberry_Gene")

    # NOTE, this is the merger of the blueberry genes, their module IDs, and
    # the multiple log2FC files received from Melanie.
    all_merged_float = pd.merge(
        log_merged_float_format,
        blueberry_genes_and_module_colors_table,
        on="Blueberry_Gene",
        how="outer",
    )

    # NB add the number of genes per module to the main table
    all_merged_float = all_merged_float.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )

    # NOTE
    # Convert positive values in Melanie's log2FC data to a new col with
    # string 'Up' and take the negative values and convert them to a new col
    # with string 'Down'. This will help in performing a value_counts() function
    # later. Take the old column of float values and convert to the string
    # 'Total_Regulated'.

    # TODO put this in a function or something
    for column in expression_contexts:
        all_merged_float.loc[
            all_merged_float[column] < 0, column + "_Direction"
        ] = "Down"
        all_merged_float.loc[all_merged_float[column] > 0, column + "_Direction"] = "Up"
        all_merged_float.loc[
            all_merged_float[column].notna(), column
        ] = "Total_Regulated"

    # NB perform the value counts
    counts_per_module = count_direction_of_regulation_per_module_log2fc(
        all_merged_float
    )

    # Get the gene counts per module information again because it was lost in
    # the previous function.
    counts_per_module = counts_per_module.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )
    counts_per_module.set_index("Module_Color", inplace=True)

    # Divide to get the percentage
    percentages = counts_per_module.apply(
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
    percentages.fillna(np.nan, inplace=True)
    percentages.set_index("Module_Color", inplace=True)
    top_n_modules_and_percent_per_context(percentages, output_dir)
    percentages.to_csv(
        os.path.join(output_dir, "Log2FC_percentages_in_modules.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

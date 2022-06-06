#!/usr/bin/env python3

"""
Filter Log2FC gene representation data
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os
from functools import reduce

from src.modules.filter_modules import read_synteny_homology_table
from src.modules.filter_modules import read_gene_modules_table


def read_log_2fc_dir(input_directory):
    # TODO no docstring and shitty documentation here
    only_log2fc_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("_FC.tsv" in f)
    ]
    if len(only_log2fc_files) != 14:
        raise ValueError

    return only_log2fc_files


def read_all_log_2fc_files(list_log_2fc_files):
    # TODO no docstring and shitty documentation here
    list_of_log_2fc_files_read = [
        pd.read_csv(x, sep="\t", header="infer").drop(
            columns=["Arabidopsis_Gene", "E_Value", "Point_of_Origin"]
        )
        for x in list_log_2fc_files
    ]
    return list_of_log_2fc_files_read


def count_direction_of_regulation_per_module(all_merged):
    # TODO no docstring and shitty documentation here
    unwanted_columns = ["Blueberry_Gene", "Module_Color"]
    columns_of_direction = all_merged.columns.to_list()
    columns_of_direction = [
        column_name
        for column_name in columns_of_direction
        if column_name not in unwanted_columns
    ]

    counts_per_module_comparison = []
    for direction_column in columns_of_direction:
        x = (
            all_merged.groupby("Module_Color")[direction_column]
            .value_counts()
            .unstack()
            .drop(columns="NA")
        )
        comparison_name = x.columns.name
        x = x.add_prefix(comparison_name + "_")
        x = x.rename_axis(None, axis=1)
        counts_per_module_comparison.append(x)
    y = reduce(
        lambda left, right: pd.merge(left, right, on="Module_Color", how="outer"),
        counts_per_module_comparison,
    )
    return y


def count_direction_of_regulation(log_2fc_file):
    value_column = [col for col in log_2fc_file.columns if "_FC" in col]
    if len(value_column) > 1:
        raise ValueError
    value_column = value_column[0]
    log_2fc_file.loc[
        log_2fc_file[value_column] > 0, [value_column + "_Direction"]
    ] = "Up"
    log_2fc_file.loc[
        log_2fc_file[value_column] < 0, [value_column + "_Direction"]
    ] = "Down"
    log_2fc_file.drop(columns=[value_column], inplace=True)  # drop the float
    # val column now that we have gotten the new column
    return log_2fc_file


def calc_gene_counts_per_module(blueberry_genes_and_module_colors_table):
    new_table = blueberry_genes_and_module_colors_table.groupby("Module_Color").count()
    new_table.rename(
        columns={"Blueberry_Gene": "Number_of_Blueberry_Genes_in_Module"}, inplace=True
    )
    return new_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # NOTE files will have to be iterated over
    parser.add_argument(
        "log_2fc_dir", type=str, help="received from Melanie on Google Drive"
    )

    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )

    parser.add_argument(
        "synteny_homology_table",
        type=str,
        help="file representing the blueberry genes and their Arabidopsis orthologs",
    )
    args = parser.parse_args()
    args.log_2fc_dir = os.path.abspath(args.log_2fc_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.synteny_homology_table = os.path.abspath(args.synteny_homology_table)
    output_path = os.path.abspath(args.log_2fc_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read the datasets
    synteny_homology_table = read_synteny_homology_table(args.synteny_homology_table)
    blueberry_genes_and_module_colors_table = read_gene_modules_table(
        args.genes_and_module_colors
    )

    gene_counts_per_module = calc_gene_counts_per_module(
        blueberry_genes_and_module_colors_table
    )
    gene_counts_per_module.reset_index(inplace=True)

    list_log_2fc_files = read_log_2fc_dir(args.log_2fc_dir)

    pandas_log_2fc_files = read_all_log_2fc_files(list_log_2fc_files)
    pandas_log_2fc_files_w_status = [
        count_direction_of_regulation(x) for x in pandas_log_2fc_files
    ]

    log_merged = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        pandas_log_2fc_files,
    )
    log_merged.fillna("NA", inplace=True)

    all_merged = pd.merge(
        log_merged,
        blueberry_genes_and_module_colors_table,
        on="Blueberry_Gene",
        how="outer",
    )
    final_counts = count_direction_of_regulation_per_module(all_merged)
    # print(gene_counts_per_module)

    x = pd.merge(final_counts, gene_counts_per_module, on="Module_Color", how="outer")
    x.set_index("Module_Color", inplace=True)

    percentages = x.apply(
        lambda row: row / row["Number_of_Blueberry_Genes_in_Module"], axis=1
    )  # NB calculate percent
    percentages.drop(
        columns=["Number_of_Blueberry_Genes_in_Module"], inplace=True
    )  # NB remove column
    percentages = percentages.add_prefix("Percent_")  # NB add prefix to column
    # names

    x = x.add_prefix("Number_Genes_")  # NB add prefix to column names for the
    # gene counts

    # NOTE at some point I need to fill nas for missing values since not all
    # genes will be 2log fc expression in a given context.

    # NOTE merge all for final result
    # TODO get arg parse results folder
    final = pd.merge(percentages, x, on="Module_Color", how="outer")
    final.to_csv("Final.tsv", sep="\t", header=True)

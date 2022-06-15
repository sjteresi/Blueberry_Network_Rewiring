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
    """
    Reads an input directory and returns a list of file paths that have the
    magic substring '_FC.tsv' in it.

    Args:
        input_directory (str): Path to directory that contains files with our
        magic substring.

    Returns:
        only_log2fc_files (list): List of file paths
    """
    only_log2fc_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("_FC.tsv" in f)
    ]  # MAGIC substring for filename recognition
    if len(only_log2fc_files) != 14:
        raise ValueError
    return only_log2fc_files


def read_all_log_2fc_files(list_log_2fc_files):
    """
    Reads a list of file paths and reads them as pandas files. Creates a list
    of pandas dataframes. Magic column names are read from the file

    Args:
        list_log_2fc_files (str): Output from read_log_2fc_dir function. List
        of file paths.

    Returns:
        list_of_log_2fc_files_read (list): List of pandas dataframes
    """
    list_of_log_2fc_files_read = [
        pd.read_csv(x, sep="\t", header="infer").drop(
            columns=["Arabidopsis_Gene", "E_Value", "Point_of_Origin"]
        )
        for x in list_log_2fc_files
    ]  # MAGIC column names for file column recognition
    return list_of_log_2fc_files_read


def count_direction_of_regulation_per_module(all_merged):
    """
    Converts the Log2FC columns of all_merged into counts per module. Index is
    now the module ID not the gene.
    """
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

    # Perform a pandas merge on a list
    y = reduce(
        lambda left, right: pd.merge(left, right, on="Module_Color", how="outer"),
        counts_per_module_comparison,
    )
    return y


def count_direction_of_regulation(log_2fc_file):
    """
    Takes a Log2FC file, (received from Melanie and read in as a pandas
    DataFrame, and determines the direction of regulation. It converts negative
    float values to 'Down' and positive float values to 'Up'.

    Args:
        log_2fc_file (pandas.DataFrame): An individual output of
        read_all_log_2fc_files

    Returns:
        log_2fc_file (pandas.DataFrame)
    """
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
    # NOTE, at this point the WGCNA data is fully read and we have the number
    # of genes per module.

    list_log_2fc_files = read_log_2fc_dir(args.log_2fc_dir)
    pandas_log_2fc_files = read_all_log_2fc_files(list_log_2fc_files)
    pandas_log_2fc_files_w_status = [
        count_direction_of_regulation(x) for x in pandas_log_2fc_files
    ]

    # NOTE, at this point I also have filtered Log2FC files from Melanie that
    # have a string for the direction (up, down)

    # NB do a pandas merge over a list, merge all the log2fc files (which have
    # different columns), but a common blueberry gene column
    log_merged = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        pandas_log_2fc_files,
    )
    log_merged.fillna("NA", inplace=True)  # NB this should fill in the
    # locations where a gene wasn't present in one comparison, but is present
    # in another. This file is the direction files all merged together

    # NOTE, now merge with the filtered WGCNA data that has been transformed to
    # have gene counts.
    all_merged = pd.merge(
        log_merged,
        blueberry_genes_and_module_colors_table,
        on="Blueberry_Gene",
        how="outer",
    )

    # NB there are some genes in the the blueberry modules dataset that were
    # never recovered in the log2fc dataset, fill those with NA. Could change
    # the 'how' keyword argument to 'left' to not add the extra genes.
    all_merged.fillna("NA", inplace=True)

    # NB key step, convert all_merged's down and up values into counts per
    # module
    direction_counts_per_context = count_direction_of_regulation_per_module(all_merged)

    direction_counts_per_context_and_gene_counts_per_module = pd.merge(
        direction_counts_per_context,
        gene_counts_per_module,
        on="Module_Color",
        how="outer",
    )
    direction_counts_per_context_and_gene_counts_per_module.set_index(
        "Module_Color", inplace=True
    )

    percentages = direction_counts_per_context_and_gene_counts_per_module.apply(
        lambda row: row / row["Number_of_Blueberry_Genes_in_Module"], axis=1
    )  # NB calculate percent, make it a new pandas dataframe
    percentages.drop(
        columns=["Number_of_Blueberry_Genes_in_Module"], inplace=True
    )  # NB remove column
    percentages = percentages.add_prefix("Percent_")  # NB add prefix to column
    # names

    final = direction_counts_per_context_and_gene_counts_per_module.add_prefix(
        "Number_Genes_"
    )  # NB add prefix to column names for the
    # gene counts

    # NOTE merge all for final result
    final = pd.merge(percentages, final, on="Module_Color", how="outer")
    final.to_csv(
        os.path.join(output_path, "Melanie_Log_2FC_Filtered.tsv"), sep="\t", header=True
    )

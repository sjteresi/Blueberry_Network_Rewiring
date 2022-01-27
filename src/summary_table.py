#!/usr/bin/env python3

"""
Create the union of dataframes of differentially expressed genes, Arabidopsis
ortholog, Arabidopsis GO term, and blueberry gene network module identity.
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


def read_TPM_expression_table(input_file):
    """
    Read a TPM expression table from file using Pandas
    """
    data = pd.read_csv(input_file, sep="\t", header="infer")
    data.rename(columns={"Gene_Name": "Blueberry_Gene"}, inplace=True)
    return data


def read_diff_ex_direction_directory(input_directory):
    """
    Parse over a directory and get the files that do not have the suffix
    "_Summary", the files I want are result files from the DEG portion of the
    project, they specify the up or down regulation of a gene in a context.

    The code then gets the comparison ID from the filename and calls a function
    to read it as a Pandas dataframe, returning a list of Pandas dataframes
    """
    only_direction_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("_Summary" not in f)
    ]

    # MAGIC get filename without extension
    comparison_groupings = [
        os.path.basename(os.path.splitext(f)[0]).rstrip("_Direction")
        for f in only_direction_files
    ]

    # NB sort alphabetically so we can zip and then pair up
    only_direction_files.sort()
    comparison_groupings.sort()

    table_list = []
    for direction_data, comparison_group in zip(
        only_direction_files, comparison_groupings
    ):
        table = read_diff_ex_direction_file(direction_data, comparison_group)
        table_list.append(table)
    return table_list


def read_diff_ex_direction_file(input_file, comparison_group):
    """
    Read a file of genes that are differentially expressed as a Pandas
    dataframe.

    Args:
        input_file (str): Path to file
        comparison_group (str): Identity of the RNA-seq libraries being
        compared

    Returns:
        pandas dataframe
    """
    data = pd.read_csv(input_file, sep="\t", header="infer")
    data.rename(
        columns={
            "Gene_Name": "Blueberry_Gene",
            "Direction_Differentially_Regulated": comparison_group,
        },
        inplace=True,
    )
    return data


def read_arabidopsis_go_table(input_file):
    """
    Read a file of Arabidopsis genes and their GO terms. Reformat the GO terms
    to be more legible.

    Args:
        input_file (str): Path to file

    Returns:
        pandas dataframe
    """
    data = pd.read_csv(
        input_file, sep="\t", header=None, names=["Arabidopsis_Gene", "GO_Terms"]
    )
    data["GO_Terms"] = [x.replace(" ", "").split(",") for x in data["GO_Terms"]]
    return data


def merge_dataframes(
    synteny_homology_data,
    arabidopsis_go_data,
    blueberry_genes_and_module_colors,
    blueberry_tpm,
    diff_ex_panda_list,
):
    """

    Args:
        synteny_homology_data (pandas.core.frame.DataFrame): Dataframe of
            blueberry genes and their Arabidopsis counterparts
        arabidopsis_go_data (pandas.core.frame.DataFrame): Dataframe of
            Arabidopsis genes and their associated GO terms
        blueberry_genes_and_module_colors (pandas.core.frame.DataFrame):
            Dataframe of blueberry genes and their module IDs (strings)
        blueberry_tpm (pandas.core.frame.DataFrame): Dataframe of blueberry
            genes and their gene expression floats
        diff_ex_panda_list (): Dataframe of blueberry genes and the direction
            (int) they are differentially expressed in each context

    Returns:
        complete (pandas.core.frame.DataFrame): Dataframe that contains all of
        the information of the dataframes above, save for genes in the GO
        dataset that were not present in the Arabidopsis ortholog set.
    """

    # NOTE merge blueberry gene / module color with blueberry gene / diff ex
    # status
    all_diff_ex_merged = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        diff_ex_panda_list,
    )
    merged = pd.merge(
        all_diff_ex_merged,
        blueberry_genes_and_module_colors,
        on="Blueberry_Gene",
        how="outer",
    )

    # NOTE merge in synteny homology data to get Arabidopsis gene, some
    # Arabidopsis genes will repeat, and there are some blueberry genes from
    # scaffolds that aren't the main chromosomes. They will have NAs in the
    # diff ex and module set.
    # This has all of the blueberry genes, no duplicates
    # SQL full outer join
    merged = pd.merge(merged, synteny_homology_data, on="Blueberry_Gene", how="outer")

    # NOTE merge in blueberry TPM data
    merged = pd.merge(merged, blueberry_tpm, on="Blueberry_Gene", how="outer")

    # NOTE get the GO data for an Arabidopsis gene only if it is already in the
    # data table, SQL left outer join
    complete = pd.merge(merged, arabidopsis_go_data, on="Arabidopsis_Gene", how="left")

    return complete


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

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

    # NOTE 'direction' files will have to be iterated over
    parser.add_argument(
        "diff_ex_dir",
        type=str,
        help="""directory containing the differential expressionfiles, files
        are divided by comparison ID and representup/down differentially
        regulated""",
    )

    parser.add_argument(
        "arabidopsis_go_terms",
        type=str,
        help="file representing an Arabidopsis gene and its GO term list",
    )

    parser.add_argument(
        "blueberry_tpm",
        type=str,
        help="""file representing a
                        blueberry gene and its TPM expression level on a
                        library by library basis""",
    )

    parser.add_argument(
        "output_filename",
        type=str,
        help="string representing the name for the output (summary table)",
    )

    parser.add_argument("output_dir", type=str, help="parent path to output results")

    args = parser.parse_args()
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.synteny_homology_table = os.path.abspath(args.synteny_homology_table)
    args.diff_ex_dir = os.path.abspath(args.diff_ex_dir)
    args.arabidopsis_go_terms = os.path.abspath(args.arabidopsis_go_terms)
    args.blueberry_tpm = os.path.abspath(args.blueberry_tpm)
    args.output_filename = os.path.abspath(args.output_filename)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read the datasets
    synteny_homology_data = read_synteny_homology_table(args.synteny_homology_table)
    arabidopsis_go_data = read_arabidopsis_go_table(args.arabidopsis_go_terms)
    blueberry_genes_and_module_colors = read_gene_modules_table(
        args.genes_and_module_colors
    )
    blueberry_tpm = read_TPM_expression_table(args.blueberry_tpm)
    diff_ex_panda_list = read_diff_ex_direction_directory(args.diff_ex_dir)

    # Merge the datasets
    complete = merge_dataframes(
        synteny_homology_data,
        arabidopsis_go_data,
        blueberry_genes_and_module_colors,
        blueberry_tpm,
        diff_ex_panda_list,
    )

    # Fill NA values with a string, because default pandas behavior is to just
    # have an empty cell
    complete.fillna("NA", inplace=True)

    # Save the merged dataset
    complete.to_csv(
        os.path.join(args.output_dir, args.output_filename),
        sep="\t",
        header=True,
        index=False,
    )

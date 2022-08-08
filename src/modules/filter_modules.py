#!/usr/bin/env python3

# TODO edit description to remove references to differential expression set and
# describe that this step is for TOPGO. Describe better that we are only
# dropping dupliccates on a module-by-module basis
"""
Execution file

Find union of differential expression / orthology set with the WGCNA output of
genes assigned to modules.
"""

__author__ = "Scott Teresi"

import argparse
import logging
import coloredlogs
import os
import pandas as pd

from src.read_tables_and_dir import read_synteny_homology_table, read_gene_modules_table


def save_color_groups(color, dataframe, output_dir):
    dataframe.drop(columns=["Module_Color", "Blueberry_Gene"], inplace=True)
    dataframe.to_csv(
        os.path.join(output_dir, str(color + "_module_AT_genes.tsv")),
        sep="\t",
        header=False,
        index=False,
    )


def save_missing_color_groups(color, dataframe, output_dir):
    dataframe.drop(columns=["Module_Color"], inplace=True)
    dataframe.to_csv(
        os.path.join(output_dir, str(color + "_module_BB_genes.tsv")),
        sep="\t",
        header=False,
        index=False,
    )


def no_match_arabidopsis(gene_modules, gene_pairs):
    """
    Identify genes that possess a module identity, that DO NOT have an
    Arabidopsis ortholog

    Args:
        gene_modules ():
        gene_pairs ():


    Returns:
        missing_genes (pandas dataframe): A pandas frame that is similar in
        structure to gene_modules, it has a Gene_Name and Module_Color column.
        And the entries in this dataframe are the genes that did not have a
        correpsonding pair in gene_pairs. So in effect we return a dataframe of
        blueberry genes in modules that do not have a corresponding AT gene.
    """
    missing_genes = gene_modules[
        ~gene_modules["Blueberry_Gene"].isin(gene_pairs.Blueberry_Gene.tolist())
    ]
    return missing_genes


def process(genes_w_module_groups, gene_pairs, output_dir):
    """

    Args:
        genes_w_module_groups (str): Path to file
        gene_pairs (str): Path to file
        output_dir (str): Path to directory
    """
    gene_modules_table = read_gene_modules_table(genes_w_module_groups)
    synteny_homology_table = read_synteny_homology_table(gene_pairs)

    blueberry_genes_w_no_match = no_match_arabidopsis(
        gene_modules_table, synteny_homology_table
    )

    # MAGIC column name
    grouped_modules = blueberry_genes_w_no_match.groupby(["Module_Color"])
    for color, frame in grouped_modules:
        output_location = os.path.join(output_dir, "missing_blueberry_modules")
        if not os.path.exists(output_location):
            os.makedirs(output_location)
        save_missing_color_groups(color, frame, output_location)

    # Now get a master table of all non-missing Arabidopsis-Blueberry pairs
    modules_w_AT = pd.merge(
        gene_modules_table, synteny_homology_table, on="Blueberry_Gene"
    )
    modules_w_AT.drop(columns=["Point_of_Origin"], inplace=True)

    # MAGIC column name
    grouped_modules = modules_w_AT.groupby(["Module_Color"])

    for color, frame in grouped_modules:
        # NOTE
        # We sort by E-value because instead of dropping all duplicates, we
        # keep the first occurrence of a duplicate. Values have been sorted by
        # E-Value such that we retain the best (lowest) E-Value score.
        frame = frame.sort_values(by=["E_Value"])
        frame = frame.drop_duplicates(subset=["Arabidopsis_Gene"], keep="first")
        frame.drop(columns=["E_Value"], inplace=True)

        output_location = os.path.join(output_dir, "modulecolors_AT")
        if not os.path.exists(output_location):
            os.makedirs(output_location)
        save_color_groups(color, frame, output_location)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="generate union of data")
    path_main = os.path.abspath(__file__)

    parser.add_argument(
        "genes_w_module_groups",
        type=str,
        help="parent path of wgcna output file containing genes and module groupings",
    )

    parser.add_argument(
        "gene_pairs",
        type=str,
        help="path to synteny/homology file containing AT BB pairs",
    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory for data output",
    )

    parser.add_argument(
        "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    args.genes_w_module_groups = os.path.abspath(args.genes_w_module_groups)
    args.gene_pairs = os.path.abspath(args.gene_pairs)
    args.output_dir = os.path.abspath(args.output_dir)

    process(
        args.genes_w_module_groups,
        args.gene_pairs,
        args.output_dir,
    )

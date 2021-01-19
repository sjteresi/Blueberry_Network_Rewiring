#!/usr/bin/env python3

"""
Master controller file. Tie together and invoke all the scripts.

Find union of differential expression / orthology set with the WGCNA output of
genes assigned to modules.
"""

__author__ = "Scott Teresi"

import argparse
import logging
import coloredlogs
import os
import pandas as pd


from scripts.modules import raise_if_no_file, raise_if_no_dir
from scripts.modules.import_modules import Gene_Modules
from scripts.modules.import_union_diffex_orth import DiffEx_Orth
from scripts.modules.import_pairs import Gene_Pairs


def validate_args(args, logger):
    """Raise if an input argument is invalid."""
    raise_if_no_file(
        args.genes_w_module_groups,
        logger=logger,
        msg_fmt="arg 'genes_w_module_groups' not a file: %s",
    )

    raise_if_no_file(
        args.gene_pairs, logger=logger, msg_fmt="arg 'gene_pairs' not a file: %s",
    )

    raise_if_no_dir(
        args.diffex_orth_dir,
        logger=logger,
        msg_fmt="arg 'diffex_orth_dir' not a dir: %s",
    )

    raise_if_no_dir(
        args.output_dir, logger=logger, msg_fmt="arg 'output_dir' not a dir: %s"
    )


def replace_wgcna_genes_w_AT(gene_modules_instance, gene_pairs_instance):
    """

    """
    merged_dataframe = pd.merge(
        gene_modules_instance.dataframe,
        gene_pairs_instance.dataframe,
        left_on="Gene_Names",
        right_on="Blueberry_Gene",
    )
    merged_dataframe.drop(
        columns=["E_Value", "Point_of_Origin", "Gene_Names"], inplace=True
    )
    return merged_dataframe


def subset_module_colors(modules_w_AT):
    """

    """
    grouped_colors = modules_w_AT.groupby(["moduleColor"])
    return grouped_colors


def save_color_groups(color, dataframe, output_dir):
    """

    """
    dataframe.drop(columns=["moduleColor", "Blueberry_Gene"], inplace=True)
    dataframe.to_csv(
        os.path.join(output_dir, str(color + "_module_AT_genes.tsv")),
        sep="\t",
        header=False,
        index=False,
    )


def save_missing_color_groups(color, dataframe, output_dir):
    """

    """
    dataframe.drop(columns=["moduleColor"], inplace=True)
    dataframe.to_csv(
        os.path.join(output_dir, str(color + "_module_BB_genes.tsv")),
        sep="\t",
        header=False,
        index=False,
    )


def drop_duplicate_AT_genes(modules_w_AT):
    """
    Removes all duplicate values for a given row, by default keeps 1 entry,
    the first occurrence of that value.
    """
    modules_w_AT.drop_duplicates(subset=["Arabidopsis_Gene"], inplace=True)
    return modules_w_AT


def no_match_blueberry(gene_modules, gene_pairs):
    """
    Args:
        gene_modules ():
        gene_pairs ():


    Returns:
        missing_genes (pandas dataframe): A pandas frame that is similar in
        structure to gene_modules, it has a Gene_Name and moduleColor column.
        And the entries in this dataframe are the genes that did not have a
        correpsonding pair in gene_pairs. So in effect we return a dataframe of
        blueberry genes in modules that do not have a corresponding AT gene.
    """
    missing_genes = gene_modules.dataframe[
        ~gene_modules.dataframe["Gene_Names"].isin(gene_pairs.blueberry_gene_list)
    ]
    return missing_genes


def process(genes_w_module_groups, gene_pairs, diffex_orth_dir, output_dir):
    """

    Args:
        genes_w_module_groups (str): Path to file
        gene_pairs (str): Path to file
        diffex_orth_dir (str): Path to directory
        output_dir (str): Path to directory
    """
    gene_modules_instance = Gene_Modules(genes_w_module_groups)
    gene_pairs_instance = Gene_Pairs(gene_pairs)

    blueberry_genes_no_match = no_match_blueberry(
        gene_modules_instance, gene_pairs_instance
    )
    grouped_modules = subset_module_colors(blueberry_genes_no_match)
    for color, frame in grouped_modules:
        frame_copy = frame.copy(deep=True)
        save_missing_color_groups(
            color,
            frame_copy,
            os.path.join(output_dir, "WGCNA_Data", "missing_blueberry_modules"),
        )
    modules_w_AT = replace_wgcna_genes_w_AT(gene_modules_instance, gene_pairs_instance)

    grouped_modules = subset_module_colors(modules_w_AT)
    for color, frame in grouped_modules:
        frame = drop_duplicate_AT_genes(frame)  # comment this if duplicate
        # genes are wanted in the module
        frame_copy = frame.copy(deep=True)
        save_color_groups(
            color, frame_copy, os.path.join(output_dir, "WGCNA_Data", "modulecolors_AT")
        )


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
        "diffex_orth_dir",
        type=str,
        help="parent directory containing directories for each library comparison",
    )

    parser.add_argument(
        "error_correction",
        type=str,
        help="what type of error correction was used (BF, FDR)",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        default=os.path.join(path_main, "../../../../", "Blueberry_Data/"),
        help="parent directory for data output",
    )

    parser.add_argument(
        "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # TODO verify if they are files or directories

    args.genes_w_module_groups = os.path.abspath(args.genes_w_module_groups)
    args.gene_pairs = os.path.abspath(args.gene_pairs)
    args.diffex_orth_dir = os.path.abspath(args.diffex_orth_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    validate_args(args, logger)

    process(
        args.genes_w_module_groups,
        args.gene_pairs,
        args.diffex_orth_dir,
        args.output_dir,
    )

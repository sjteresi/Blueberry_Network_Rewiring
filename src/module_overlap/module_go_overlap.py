#!/usr/bin/env python3

"""
Filter GO data. Determine overlap with modules from WGCNA.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

from src.modules.filter_modules import read_gene_modules_table
from src.module_overlap.module_log2fc_overlap import (
    calc_gene_counts_per_module,
)


def read_GO_dir(input_directory):
    """
    Reads an input directory and returns a list of file paths that contain
    the magic substring 'Summary'

    Args:
        input_directory (str): Path to directory that contains files with our
        magic substring. This directory is parsed to find filenames that match
        the substring.

    Returns:
        go_files (list): List of absolute file paths
    """
    go_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("Summary" in f)
    ]  # MAGIC substring for filename recognition

    return go_files


def read_GO_file(go_files):
    """
    Reads a list of file paths and reads them as pandas files. Creates a list
    of pandas dataframes. Magic column names are read from the file

    Args:
        go_files (str): Output from read_log_2fc_dir function. List
        of file paths.

    Returns:
    """
    my_dict = {}  # TODO refactor, because I am basically
    for go_file_path in go_files:
        individual_go_panda = pd.read_csv(go_file_path, sep="\t", header="infer")
        individual_go_panda.rename(
            columns={
                "GO.ID": "GO_ID",
                "Term": "GO_TERM",
            },
            inplace=True,
        )  # MAGIC column name
        individual_go_panda.drop(
            columns=[
                "Annotated",
                "Significant",
                "Expected",
                "classicFisher",
                "Overrepresented",
            ],
            inplace=True,
        )

        if len(os.path.splitext(os.path.basename(go_file_path))) > 2:
            # MAGIC remove filename extension, does not work if multiple periods
            # are in filename.
            raise ValueError(
                """There were multiple periods in your file path,
                             this will cause an error."""
            )
        # MAGIC to get filename
        analysis_context = os.path.splitext(os.path.basename(go_file_path))[0]

        # MAGIC string split on underscores in filename, there is a common
        # structure to the filename that allows me to do this
        module_name = analysis_context.split("_")[0]
        analysis_type = analysis_context.split("_")[4]

        individual_go_panda["GO_Term_Type"] = analysis_type
        individual_go_panda[module_name] = "Recovered"

        # TODO refactor
        if module_name not in my_dict:
            my_dict[module_name] = []
        my_dict[module_name].append(individual_go_panda)

    return my_dict


def read_interesting_GO_terms(input_path):
    colnames = ["GO_ID", "GO_Term"]
    interesting_go_terms = pd.read_csv(
        input_path, names=colnames, header=None, sep="\t"
    )
    return interesting_go_terms


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("go_dir", type=str, help="output from topGO")
    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )
    parser.add_argument(
        "interesting_go_terms",
        type=str,
        help="file representing the interesting GO terms received from Melanie",
    )
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.go_dir = os.path.abspath(args.go_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.interesting_go_terms = os.path.abspath(args.interesting_go_terms)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    go_file_paths = read_GO_dir(args.go_dir)
    the_dictionary = read_GO_file(go_file_paths)

    # Need to concat the files if they have the same columns
    all_go_list = [pd.concat(val) for key, val in the_dictionary.items()]

    # Perform a pandas merge on a list
    all_go_merged = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=["GO_ID", "GO_TERM", "GO_Term_Type"],
            how="outer",
        ),
        all_go_list,
    )

    # NOTE
    # Remove rows where the entire row is NA
    # Remove columns (modules) where the entire module is NA
    all_go_merged.dropna(axis=0, how="all", inplace=True)
    all_go_merged.dropna(axis=1, how="all", inplace=True)
    all_go_merged.fillna("Not_Recovered", inplace=True)

    interesting_go_terms = read_interesting_GO_terms(args.interesting_go_terms)

    interesting_go_term_info = all_go_merged[
        all_go_merged["GO_ID"].isin(interesting_go_terms["GO_ID"])
    ]

    interesting_go_terms = interesting_go_terms.sort_values(
        by=["GO_ID"], ascending=True
    )

    interesting_go_term_info.to_csv(
        os.path.join(args.output_dir, "Interesting_GO_Term_Module_Representation.tsv"),
        sep="\t",
        header=True,
        index=False,
    )

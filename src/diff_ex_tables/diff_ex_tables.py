#!/usr/bin/env python3
# NOTE
# Candidate for deletion, constructed table by hand.

"""
TODO fill out
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
from collections import namedtuple
import logging
import os
import pandas as pd

DifferentialCounts = namedtuple(
    "DifferentialCounts", ["Identity", "Comparison_Order", "Up", "NotSig", "Down"]
)


def load_differential_expression_data(diff_ex_dir):
    """
    Load the data

    Args:
        diff_ex_dir (str): Path to folder containing tsv files of gene
        identities
    """
    list_of_summary_files = []
    identities = []
    for individual_file in os.scandir(diff_ex_dir):
        if individual_file.path.endswith("Summary.txt"):
            list_of_summary_files.append(os.path.join(diff_ex_dir, individual_file))

    list_of_up_down_gene_counts = []
    for individual_summary_file in list_of_summary_files:
        summary_info = pd.read_csv(
            individual_summary_file,
            sep="\t",
            header=None,
            names=["Comparison_Order", "Count"],
        )
        identity = os.path.basename(individual_summary_file)
        identity = identity.strip("_Summary.txt")
        comparison_order = summary_info.iloc[0].Comparison_Order
        summary_info = summary_info.iloc[1:]
        summary_info.set_index("Comparison_Order", inplace=True)
        summary_info = summary_info.T
        summary_info.reset_index(drop=True, inplace=True)
        summary_info = summary_info.astype("int32")
        X = DifferentialCounts(
            identity,
            comparison_order,
            summary_info.iloc[0]["Down"],
            summary_info.iloc[0]["NotSig"],
            summary_info.iloc[0]["Up"],
        )


def create_table():
    raise NotImplementedError


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # parser.add_argument("", type=str, help="")
    parser.add_argument("diff_ex_dir", type=str, help="TODO")
    parser.add_argument("output_dir", type=str, help="parent path to output results")
    args = parser.parse_args()
    args.diff_ex_dir = os.path.abspath(args.diff_ex_dir)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    load_differential_expression_data(args.diff_ex_dir)

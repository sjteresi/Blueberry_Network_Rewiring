#!/usr/bin/env python3

"""
Calculates the proportions of each gene in the genome belonging to each
identification type (Synteny, BLAST, Not found)
"""

__author__ = "Alder Fulton, Scott Teresi"

import argparse
import os
import pandas as pd
import numpy as np


def process(gene_annotation, orthology_file, efilter):
    """
    Args:
       gene_annotation (str): path to the csv containing the data for different genes.
       orthology_file: path to the orthology file containg the different identification types.

    Returns:
        None. Prints the proportions of each gene.

    """
    # Import gene_annotation
    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stope",
        "Score",
        "Strand",
        "Frame",
        "Fullname",
    ]

    col_to_use = [
        "Feature",
        "Fullname",
    ]

    genes_df = pd.read_csv(
        gene_annotation,
        delimiter="\t",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype="str",
    )
    genes_df = genes_df[genes_df["Feature"] == "gene"]
    genes_df["Fullname"] = genes_df["Fullname"].str.split("Name=").str.get(1)
    N = genes_df["Fullname"].size

    # Import orthology_file
    orthology_df = pd.read_csv(
        orthology_file,
        delimiter="\t",
        header="infer",
        engine="python",
    )
    orthology_df.drop(columns=["Arabidopsis_Gene"], inplace=True)
    # Column name is E_Value
    orthology_df = orthology_df[orthology_df["E_Value"] <= efilter]
    S = orthology_df["Point_of_Origin"][
        orthology_df["Point_of_Origin"] == "Synteny"
    ].size
    B = orthology_df["Point_of_Origin"][orthology_df["Point_of_Origin"] == "BLAST"].size
    Q = N - (B + S)
    proportions = calc_percents(N, (S, B, Q))

    print("\nTotal Genes = " + str(N) + "\n")
    print("\nE value cutoff = " + str(efilter) + "\n")
    print(
        "Percent synteny = "
        + proportions[0]
        + "\nPercent blast = "
        + proportions[1]
        + "\nPercent neither = "
        + proportions[2]
    )


def calc_percents(whole, parts):
    """
    args:
        whole (int): the thing you're calculating the percentages of.
        parts (tuple): the parts of the whole.
    returns:
        sol (list): A list of the solutions in the same order as parts
    """
    sol = []
    for part in parts:
        sol.append(str(round(part / whole, 3) * 100) + "%")
    return sol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate gene proportions")
    path_main = os.path.abspath(
        __file__
    )  # The path that leads you to this file, not including this file.

    parser.add_argument(
        "gene_annotation",
        type=str,
        help="parent path to the gene annotation file",
    )
    parser.add_argument(
        "orthology_file",
        type=str,
        help="parent path to the synteny/BLAST file",
    )
    parser.add_argument(
        "--e_filter",
        type=float,
        help="The lower bound of the probability that the measurement was a coincidence that you want to filter out.",
        default=1e-5,
    )
    args = parser.parse_args()
    args.gene_annotation = os.path.abspath(args.gene_annotation)
    args.orthology_file = os.path.abspath(args.orthology_file)

    process(args.gene_annotation, args.orthology_file, args.e_filter)

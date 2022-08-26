#!/usr/bin/env python3

"""
Filter the master Arabidopsis TAIR protein table so that Melanie has an easier
table to work with.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os

# from functools import reduce

# from src.read_tables_and_dir import read_synteny_homology_table, read_gene_modules_table


def remove_sequences(filepath, output_dir):
    """
    Verbose way just to remove all of the protein sequence lines
    """
    half_filtered_file = os.path.join(output_dir, "Protein_and_Genes_Unfiltered.txt")
    print(half_filtered_file)
    with open(filepath, "r") as fh:
        with open(half_filtered_file, "w") as new_file:
            for curline in fh:
                if not curline.startswith(">"):
                    continue
                else:
                    curline = curline.lstrip(">")
                    new_file.write(curline)
    return half_filtered_file


def read_proteins(filepath, output_dir):
    col_names = ["Arabidopsis_Gene", "Protein_ID", "Protein_Name", "Location"]
    data = pd.read_csv(filepath, sep="|", skipinitialspace=True, names=col_names)
    data.drop(columns="Location", inplace=True)

    # Subsitute and remove various weird strings in data set to make legible
    data["Protein_ID"] = data["Protein_ID"].str.replace(" ", "")
    data["Protein_ID"] = data["Protein_ID"].str.lstrip("Symbols:")
    data.loc[data["Protein_ID"] == "", "Protein_ID"] = "No_ID"

    # data.loc[data["Protein_ID"] == "no symbol available ", "Protein_ID"] = "NA"
    # data.loc[data["Protein_Name"] == "no full name available ", "Protein_Name"] = "NA"

    # NB split the Arabidopsis gene names by the dot notation and only keep the
    # first part.
    data["Arabidopsis_Gene"] = data["Arabidopsis_Gene"].str.split(".").str[0]

    # NB drop duplicate rows if ALL columns are the same
    data = data.drop_duplicates()
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "arabidopsis_protein_file",
        type=str,
        help="""master file derived from TAIR, contains gene protein ID and
        function info""",
    )
    parser.add_argument(
        "output_dir", type=str, help="string representing the name for the output dir"
    )
    args = parser.parse_args()
    args.arabidopsis_protein_file = os.path.abspath(args.arabidopsis_protein_file)
    args.output_dir = os.path.abspath(args.output_dir)

    half_filtered_file = remove_sequences(
        args.arabidopsis_protein_file, args.output_dir
    )
    data = read_proteins(half_filtered_file, args.output_dir)
    data.to_csv(
        os.path.join(args.output_dir, "Filtered_Arabidopsis_Protein_Info.tsv"),
        index=False,
        sep="\t",
        header=True,
    )

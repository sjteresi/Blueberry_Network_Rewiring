#!/usr/bin/env python3

"""
Filter the master Arabidopsis TAIR protein table so that Melanie has an easier
table to work with. Generates a table of Arabidopsis genes and their protein ID
and protein names

Writes:
    Protein_and_Genes_Unfiltered.txt (TEMP)
    Filtered_Arabidopsis_Protein_Info.tsv
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os


def remove_sequences(filepath, output_dir):
    """
    Verbose way just to remove all of the protein sequence lines, the lines
    starting with '>'

    Args:
        filepath (str)

        output_dir (str)
    """
    half_filtered_file = os.path.join(output_dir, "Protein_and_Genes_Unfiltered.txt")
    with open(filepath, "r") as fh:
        with open(half_filtered_file, "w") as new_file:
            for curline in fh:
                if not curline.startswith(">"):
                    continue
                else:
                    curline = curline.replace(">", "")
                    new_file.write(curline)
    return half_filtered_file


def filter_protein_table(filepath, output_dir):
    """
    Filter the master Arabidopsis TAIR protein/gene table with Pandas into a
    more manageable format

    Args:
        filepath (str)

        output_dir (str)
    """
    col_names = ["Arabidopsis_Gene", "Protein_ID", "Protein_Name", "Location"]
    data = pd.read_csv(filepath, sep="|", skipinitialspace=True, names=col_names)
    data.drop(columns="Location", inplace=True)

    # Subsitute and remove various weird strings in data set to make legible
    data["Protein_ID"] = data["Protein_ID"].str.replace(" ", "")
    data["Protein_ID"] = data["Protein_ID"].str.replace("Symbols:", "")
    data.loc[data["Protein_ID"] == "", "Protein_ID"] = "No_ID"

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
    data = filter_protein_table(half_filtered_file, args.output_dir)
    data.to_csv(
        os.path.join(args.output_dir, "Filtered_Arabidopsis_Protein_Info.tsv"),
        index=False,
        sep="\t",
        header=True,
    )

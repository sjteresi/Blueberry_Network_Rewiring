#!/usr/bin/env python3

"""
Test
"""

__author__ = "Scott Teresi"

import pandas as pd


def import_orthologs(syntelog_input_file):
    """

    """
    col_names = [
        "Query",
        "Subject",
        "Percent ID",
        "Alignment Length",
        "Mismatches",
        "Gap Openings",
        "Q_Start",
        "Q_Stop",
        "Subject_Start",
        "Subject_Stop",
        "E_Value",
        "Bit_Score",
    ]

    col_to_use = [
        "Query",
        "Subject",
        "E_Value",
    ]

    gene_data = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )
    gene_data.E_Value = gene_data.E_Value.astype("float64")

    # Get the correct name for the arabidopsis genes
    gene_data["Subject"] = gene_data["Subject"].str.split("\|\|").str[3]

    # Get the correct name for the blueberry genes
    gene_data["Query"] = gene_data["Query"].str.split("-mRNA-1").str[0]

    # Trim E-values less than 0.05
    gene_data = gene_data.loc[gene_data["E_Value"] < 0.05]

    # Need to take first occurrence of the gene, the one with the smallest
    # E-Value
    gene_data = gene_data.loc[gene_data.groupby("Query")["E_Value"].idxmin()]

    # Rename columns
    gene_data.rename(
        columns={"Query": "Blueberry", "Subject": "Arabidopsis"}, inplace=True
    )
    gene_data["Point_of_Origin"] = "BLAST"

    return gene_data

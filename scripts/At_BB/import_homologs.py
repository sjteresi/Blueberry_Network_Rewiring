#!/usr/bin/env python3

"""
Import the homolog data and provide a class for its access
"""

__author__ = "Scott Teresi"

import logging
import pandas as pd


def import_homologs(homolog_input_file):
    """
    Import the homologs from the raw file and manage data filtration
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
        homolog_input_file,
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

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    gene_data["Point_of_Origin"] = "BLAST"

    return gene_data

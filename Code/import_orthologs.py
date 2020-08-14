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

    Gene_Data = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )

    # Set the correct data types
    # Gene_Data.OrgA_Chromosome = Gene_Data.OrgA_Chromosome.astype(str)
    # Gene_Data.OrgB_Chromosome = Gene_Data.OrgB_Chromosome.astype(str)
    # Gene_Data.OrgA_Gene_Region = Gene_Data.OrgA_Gene_Region.astype(str)
    # Gene_Data.OrgB_Gene_Region = Gene_Data.OrgB_Gene_Region.astype(str)
    # Gene_Data.E_Value = Gene_Data.E_Value.astype("float64")
    # Gene_Data.Diagonal_Score = Gene_Data.Diagonal_Score.astype("int32")

    # Get the correct name for the arabidopsis genes
    Gene_Data["Subject"] = Gene_Data["Subject"].str.split("\|\|").str[3]

    # Get the correct name for the blueberry genes
    Gene_Data["Query"] = Gene_Data["Query"].str.split("-mRNA-1").str[0]

    # Trim E-values less than 0.05
    Gene_Data = Gene_Data.loc[Gene_Data["E_Value"] < 0.05]

    # Need to take first occurrence of the gene, the one with the smallest
    # E-Value
    Gene_Data = Gene_Data.groupby(["Query"])["E_Value"].min()

    return Gene_Data

#!/usr/bin/env python3

"""
Test
"""

__author__ = "Scott Teresi"

import pandas as pd


def import_syntelogs(syntelog_input_file):
    """

    """
    col_names = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgA_Start",
        "OrgA_Stop",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "OrgB_Start",
        "OrgB_Stop",
        "E_Value",
        "Diagonal_Score",
        "Web_Link",
    ]

    col_to_use = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "E_Value",
        "Diagonal_Score",
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
    Gene_Data.OrgA_Chromosome = Gene_Data.OrgA_Chromosome.astype(str)
    Gene_Data.OrgB_Chromosome = Gene_Data.OrgB_Chromosome.astype(str)
    Gene_Data.OrgA_Gene_Region = Gene_Data.OrgA_Gene_Region.astype(str)
    Gene_Data.OrgB_Gene_Region = Gene_Data.OrgB_Gene_Region.astype(str)
    Gene_Data.E_Value = Gene_Data.E_Value.astype("float64")
    Gene_Data.Diagonal_Score = Gene_Data.Diagonal_Score.astype("int32")

    # Get the correct name for the genes
    Gene_Data["OrgA_Gene_Region"] = (
        Gene_Data["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    Gene_Data["OrgB_Gene_Region"] = (
        Gene_Data["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    # Get the correct name for the chromosome
    Gene_Data["OrgA_Chromosome"] = Gene_Data["OrgA_Chromosome"].str.split("_").str[1]
    Gene_Data["OrgB_Chromosome"] = Gene_Data["OrgB_Chromosome"].str.split("_").str[1]

    # Trim E-values less than 0.05
    Gene_Data = Gene_Data.loc[Gene_Data["E_Value"] < 0.05]

    return Gene_Data

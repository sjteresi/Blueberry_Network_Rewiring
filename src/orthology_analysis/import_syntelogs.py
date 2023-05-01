#!/usr/bin/env python3

__author__ = "Scott Teresi"

import logging
import pandas as pd


def import_syntelogs(syntelog_input_file):
    """
    Import the syntelogs from the raw file and manage data filtration
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

    gene_data = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )

    # Set the correct data types
    gene_data.OrgA_Chromosome = gene_data.OrgA_Chromosome.astype(str)
    gene_data.OrgB_Chromosome = gene_data.OrgB_Chromosome.astype(str)
    gene_data.OrgA_Gene_Region = gene_data.OrgA_Gene_Region.astype(str)
    gene_data.OrgB_Gene_Region = gene_data.OrgB_Gene_Region.astype(str)
    gene_data.E_Value = gene_data.E_Value.astype("float64")
    gene_data.Diagonal_Score = gene_data.Diagonal_Score.astype("int32")

    # Get the correct name for the genes
    gene_data["OrgA_Gene_Region"] = (
        gene_data["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    gene_data["OrgB_Gene_Region"] = (
        gene_data["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    # Get the correct name for the blueberry genes
    gene_data["OrgB_Gene_Region"] = (
        gene_data["OrgB_Gene_Region"].str.split("-mRNA-1").str[0]
    )

    # Get the correct name for the chromosome
    gene_data["OrgA_Chromosome"] = gene_data["OrgA_Chromosome"].str.split("_").str[1]
    gene_data["OrgB_Chromosome"] = gene_data["OrgB_Chromosome"].str.split("_").str[1]

    # This step is important, it could differ if your data input is different.
    gene_data.rename(
        columns={"OrgA_Gene_Region": "Arabidopsis", "OrgB_Gene_Region": "Blueberry"},
        inplace=True,
    )
    # Trim E-values less than 0.05
    gene_data = gene_data.loc[gene_data["E_Value"] < 0.05]

    gene_data.drop(
        columns=["OrgA_Chromosome", "OrgB_Chromosome", "Diagonal_Score"],
        inplace=True,
    )

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    gene_data["Point_of_Origin"] = "Synteny"

    # Sort by Name and E-Values
    gene_data.sort_values(
        by=["Blueberry", "E_Value"], ascending=(True, True), inplace=True
    )
    # Need to take first occurrence of a duplicated gene, the one with the smallest
    # E-Value
    gene_data = gene_data.drop_duplicates(subset=["Blueberry"], keep="first")

    return gene_data

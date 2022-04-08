#!/usr/bin/env python3

__author__ = "__Alder Fulton__"

import argparse
import pandas as pd
import numpy as np
import os

# genes_colors_path = "/home/alder/Documents/research/Blueberry_Data/"
# data_folder = "/home/alder/Documents/research/Blueberry_Data/Diff_Ex/EdgeR_Output/"
# hap = "All_Hap"
# method = "Bonferroni"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adds up the differentially\
                                     expressed genes in each module for each\
                                     comparison."
    )
    parser.add_argument(
        "genes_colors_path",
        type=str,
        help="The path to the\
                        list of genes and their corresponding module colors.",
    )
    parser.add_argument(
        "data_folder_path",
        type=str,
        help="The path to the\
                        folder with all of the tsv files in it.",
    )
    parser.add_argument("hap_type", type=str, help="Single or All")
    parser.add_argument("method", type=str, help="Bonferonni or FDR")
    parser.add_argument("output_dir", type=str, help="output dir")
    args = parser.parse_args()
    args.output_dir = os.path.abspath(args.output_dir)
    args.genes_colors_path = os.path.abspath(args.genes_colors_path)
    args.data_folder_path = os.path.abspath(args.data_folder_path)

    # NB, not that good, data folder where all the diff ex data is.
    combined_path = (
        args.data_folder_path + "/" + args.hap_type + "/" + args.method + "/"
    )
    # Read the genes with module group WGCNA output
    data = pd.read_csv(args.genes_colors_path, header="infer", delimiter="\t")

    # Create the result dataframe
    result = data.groupby("moduleColor", as_index=False).count()
    result.rename(
        {"moduleColor": "Module ID", "Gene_Names": "Number of Total Genes"},
        axis=1,
        inplace=True,
    )

    for f in os.listdir(combined_path):
        # Only use tsv files
        if f[-4:] != ".tsv":
            continue

        # Reset the data for the new iteration
        genes_colors = data.copy()
        genes_colors2 = data.copy()

        # Read in the comparison data
        df = pd.read_csv(combined_path + f, delimiter="\t")

        # Replace the negatives since we only care about if it is differentially
        # expressed or not.
        genes_colors["Diff Ex"] = df["Direction_Differentially_Regulated"].replace(
            {-1: 1}
        )

        result["Diff_Ex_" + f[:-4]] = (
            genes_colors[["moduleColor", "Diff Ex"]]
            .groupby("moduleColor", as_index=False)
            .sum()["Diff Ex"]
        )

    # MAGIC filename
    result.to_csv(
        os.path.join(args.output_dir, "modules_diff_ex_table.tsv"),
        header=True,
        index=False,
        sep="\t",
    )
    # TODO why is this here?
    result["total"] = result.drop(columns=["Module ID", "Number of Total Genes"]).apply(
        np.sum, axis=1
    )

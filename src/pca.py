#!/usr/bin/env python3

"""
TODO
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import numpy as np
import pandas as pd
import os

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

from src.read_tables_and_dir import read_FPKM_or_TPM

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    # NB input args correspond to Makefile, must be ordered correctly
    parser.add_argument("tpm_file", type=str, help="tpm master file")
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.tpm_file = os.path.abspath(args.tpm_file)
    args.output_dir = os.path.abspath(args.output_dir)

    tpm_data = read_FPKM_or_TPM(args.tpm_file)
    tpm_data.set_index("Blueberry_Gene", inplace=True)
    original_col = tpm_data.columns.tolist()
    lib_c_rules = [
        "Lib_1dpo_c",
        "Lib_2dpo_c",
        "Lib_3dpo_c",
        "Lib_4dpo_c",
        "Lib_5dpo_c",
        "Lib_6dpo_c",
        "Lib_7dpo_c",
    ]
    lib_t_rules = [
        "Lib_1dpo_t",
        "Lib_2dpo_t",
        "Lib_3dpo_t",
        "Lib_4dpo_t",
        "Lib_5dpo_t",
        "Lib_6dpo_t",
        "Lib_7dpo_t",
    ]
    dra_c_rules = [
        "Dra_1dpo_c",
        "Dra_2dpo_c",
        "Dra_3dpo_c",
        "Dra_4dpo_c",
        "Dra_5dpo_c",
        "Dra_6dpo_c",
        "Dra_7dpo_c",
    ]
    dra_t_rules = [
        "Dra_1dpo_t",
        "Dra_2dpo_t",
        "Dra_3dpo_t",
        "Dra_4dpo_t",
        "Dra_5dpo_t",
        "Dra_6dpo_t",
        "Dra_7dpo_t",
    ]
    all_rules = [
        item
        for sublist in [lib_c_rules, lib_t_rules, dra_c_rules, dra_t_rules]
        for item in sublist
    ]

    for major_rule in all_rules:
        col_whitelist = [col for col in original_col if col.startswith(major_rule)]

        if major_rule.endswith("c"):
            major_rule = major_rule.replace("c", "CONTROL")
        if major_rule.endswith("t"):
            major_rule = major_rule.replace("t", "TREATMENT")

        tpm_data[major_rule + "_MEAN"] = tpm_data[col_whitelist].mean(axis=1)

    tpm_data.drop(columns=original_col, inplace=True)

    print(tpm_data)
    raise ValueError
    ############################
    tpm_data.reset_index(inplace=True)  # NOTE, this could be future problem,
    # print(tpm_data)
    tpm_data = pd.melt(tpm_data, id_vars=["Blueberry_Gene"])
    # print(tpm_data)

    # Separate features and status/target from data
    x = tpm_data.drop(columns=["Blueberry_Gene"])
    x = x.loc[:, ["value"]].values
    y = tpm_data.loc[:, ["variable"]].values

    # Start scaling the data
    scaled_tpm_data = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(scaled_tpm_data)
    principalDF = pd.DataFrame(data=principalComponents, columns=["PC_1", "PC_2"])

    tpm_data.reset_index(inplace=True)  # NOTE, this could be future problem,
    # we lost the index information when we do concat, this step is necessary
    # for the final df

    # TODO create a tag column, so I can color the dots
    final_df = pd.concat([principalDF, tpm_data[["Status"]]], axis=1)
    print(final_df)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Principal Component 1", fontsize=15)
    ax.set_ylabel("Principal Component 2", fontsize=15)
    ax.set_title("2 component PCA", fontsize=20)
    targets = final_df["Status"].unique()
    colors = ["r", "g", "b"]
    for target, color in zip(targets, colors):
        indicesToKeep = final_df["Status"] == target
        ax.scatter(
            final_df.loc[indicesToKeep, "PC_1"],
            final_df.loc[indicesToKeep, "PC_2"],
            c=color,
            s=50,
        )
    ax.legend(targets)
    ax.grid()
    plt.show()

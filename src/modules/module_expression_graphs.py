#!/usr/bin/env python3

"""
Generate linegraphs of gene expression per module
"""

__author__ = "Scott Teresi"

import argparse
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from functools import reduce

from src.read_tables_and_dir import read_gene_modules_table, read_FPKM_or_TPM
from src.module_overlap.module_log2fc_overlap import (
    read_log_2fc_dir,
    read_all_log_2fc_files,
)

graph_groups = {
    "Draper Control": [
        "Dra_1dpo_Control_Mean",
        "Dra_2dpo_Control_Mean",
        "Dra_3dpo_Control_Mean",
        "Dra_4dpo_Control_Mean",
        "Dra_5dpo_Control_Mean",
        "Dra_6dpo_Control_Mean",
        "Dra_7dpo_Control_Mean",
    ],
    "Draper Treatment": [
        "Dra_1dpo_Treatment_Mean",
        "Dra_2dpo_Treatment_Mean",
        "Dra_3dpo_Treatment_Mean",
        "Dra_4dpo_Treatment_Mean",
        "Dra_5dpo_Treatment_Mean",
        "Dra_6dpo_Treatment_Mean",
        "Dra_7dpo_Treatment_Mean",
    ],
    "Liberty Treatment": [
        "Lib_1dpo_Treatment_Mean",
        "Lib_2dpo_Treatment_Mean",
        "Lib_3dpo_Treatment_Mean",
        "Lib_4dpo_Treatment_Mean",
        "Lib_5dpo_Treatment_Mean",
        "Lib_6dpo_Treatment_Mean",
        "Lib_7dpo_Treatment_Mean",
    ],
    "Liberty Control": [
        "Lib_1dpo_Control_Mean",
        "Lib_2dpo_Control_Mean",
        "Lib_3dpo_Control_Mean",
        "Lib_4dpo_Control_Mean",
        "Lib_5dpo_Control_Mean",
        "Lib_6dpo_Control_Mean",
        "Lib_7dpo_Control_Mean",
    ],
    "All": [
        "All_1dpo_Data_Mean",
        "All_2dpo_Data_Mean",
        "All_3dpo_Data_Mean",
        "All_4dpo_Data_Mean",
        "All_5dpo_Data_Mean",
        "All_6dpo_Data_Mean",
        "All_7dpo_Data_Mean",
    ],
}


def gen_whitelist(expression_table_columns):
    """
    Generates a dictionary of expression columns in the raw TPM data and
    assigns the columns a context key. For example, Lib_6dpo_c2_S65_L007 gets
    assigned to 'Lib 6 dpo Control' and is grouped with the other RNA-seq
    replicates.
    """
    whitelist = {}
    for i in range(1, 8, 1):
        draper_control_key = "Dra_" + str(i) + "dpo_" + "Control"
        draper_treatment_key = "Dra_" + str(i) + "dpo_" + "Treatment"
        liberty_control_key = "Lib_" + str(i) + "dpo_" + "Control"
        liberty_treatment_key = "Lib_" + str(i) + "dpo_" + "Treatment"
        all_data_key = "All_" + str(i) + "dpo_" + "Data"
        draper_treatments = [
            col
            for col in expression_table_columns
            if ("Dra" in col) and (str(i) + "dpo" in col) and ("t" in col)
        ]
        draper_controls = [
            col
            for col in expression_table_columns
            if ("Dra" in col) and (str(i) + "dpo" in col) and ("c" in col)
        ]
        liberty_treatments = [
            col
            for col in expression_table_columns
            if ("Lib" in col) and (str(i) + "dpo" in col) and ("t" in col)
        ]
        liberty_controls = [
            col
            for col in expression_table_columns
            if ("Lib" in col) and (str(i) + "dpo" in col) and ("c" in col)
        ]
        all_data = [col for col in expression_table_columns if str(i) + "dpo" in col]
        whitelist[draper_control_key] = draper_controls
        whitelist[draper_treatment_key] = draper_treatments
        whitelist[liberty_control_key] = liberty_controls
        whitelist[liberty_treatment_key] = liberty_treatments
        whitelist[all_data_key] = all_data
    return whitelist


def get_context_means(whitelist, expression):
    """
    Generates a new column for the gene expression mean of each RNA-seq
    grouping
    """
    for key, val in whitelist.items():
        expression[key + "_Mean"] = expression[val].mean(axis=1)
        expression = expression.reindex(sorted(expression.columns), axis=1)
    return expression


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "master_gene_module_table",
        type=str,
        help="""parent path of wgcna output file containing
                genes and module groupings
                """,
    )
    parser.add_argument("gene_expression_table", type=str)

    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.master_gene_module_table = os.path.abspath(args.master_gene_module_table)
    args.gene_expression_table = os.path.abspath(args.gene_expression_table)
    args.output_dir = os.path.abspath(args.output_dir)

    genes_and_modules = read_gene_modules_table(args.master_gene_module_table)
    expression = read_FPKM_or_TPM(args.gene_expression_table)

    # Decipher the various RNA-seq library names to later coalesce into the
    # true RNA-seq means
    whitelist = gen_whitelist(expression.columns.to_list())
    # Mutate expression pandaframe to have the column means
    expression = get_context_means(whitelist, expression)
    for key, val in whitelist.items():
        if "All" not in key:
            expression.drop(columns=val, inplace=True)

    # Save the dataframe of expression means to disk
    expression.set_index("Blueberry_Gene").to_csv(
        os.path.join(args.output_dir, "Mean_Expression_TPM.tsv"),
        sep="\t",
        header=True,
        index=True,
    )
    # Add in columns of 'Module_Color'
    genes_modules_and_expression = pd.merge(
        expression, genes_and_modules, how="inner", on="Blueberry_Gene"
    )
    genes_modules_and_expression.set_index("Blueberry_Gene", inplace=True)

    # NOTE, MAGIC these are Module IDs that we identified a priori
    modules_of_interest = ["darkseagreen3", "plum3", "palevioletred2", "lightpink3"]

    for grouping_title, exp_columns in graph_groups.items():
        # NOTE MAGIC, split the first element of the group into pieces and tie
        # them together to get a name for the group.

        # FUTURE, melanie may want this changed later
        if "Lib" in grouping_title:
            genome_color = "red"
        elif "Dra" in grouping_title:
            genome_color = "blue"
        elif "All" in grouping_title:
            genome_color = "purple"
        else:
            raise ValueError("Group title error")

        for module in modules_of_interest:
            subsetted_by_module = genes_modules_and_expression[
                genes_modules_and_expression["Module_Color"] == module
            ].copy(deep=True)
            subsetted_by_module.drop(columns=["Module_Color"], inplace=True)
            subsetted_by_module = subsetted_by_module.filter(items=exp_columns)

            transposed = subsetted_by_module.transpose()
            transposed.sort_index(inplace=True)

            # MAGIC splitting allows us to get the day val on the 1st index
            transposed.index = [
                day_label.split("_")[1] for day_label in transposed.index.to_list()
            ]
            transposed["Mean_Gene"] = transposed.mean(axis=1)
            transposed = transposed.filter(items=["Mean_Gene"])
            log_transposed = transposed.apply(np.log10)
            plt.plot(log_transposed.index, log_transposed, color=genome_color)
            plt.xticks(log_transposed.index)

            # MAGIC ylim and yticks for visual purposes
            plt.ylim(1.3, 2.5)
            plt.yticks(np.linspace(1.3, 2.5, num=13, dtype="float32", endpoint=True))
            plt.ylabel("Mean Log10 Gene Expression (TPM)")
            title = grouping_title + "\n" + "Module ID: " + module
            plt.title(title)
            fig = plt.gcf()
            fig.set_size_inches(12, 9)
            filename = module + "_" + grouping_title.replace(" ", "_") + ".png"
            plt.savefig(
                os.path.join(
                    args.output_dir,
                    filename,
                ),
                bbox_inches="tight",
            )
            plt.clf()

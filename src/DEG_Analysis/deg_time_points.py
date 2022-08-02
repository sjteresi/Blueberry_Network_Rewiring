#!/usr/bin/env python3

"""
Create a line plot of DEGs over the time points
Create a Venn Diagram of DEGs over the time points to see overlap between
Liberty and Draper datasets
"""

__author__ = "Scott Teresi"

# General imports
import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

# Graphing functions
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted

# Import reader functions
from src.read_tables_and_dir import (
    read_DEG_dir,
    read_all_DEG_direction_files,
    read_GO_ID_w_term,
    read_synteny_homology_table,
)


def get_unique_and_non_unique_genes_per_comparison(
    all_degs_pandas, draper_col, liberty_col
):
    """
    # TODO add the arabidopsis gene column nmae

    Gets the counts of unique and non-unique DEGs per genotype comparison.
    Basically does Draper C1 and T1 vs Liberty C1 and T1. I am grouping the
    DEGs of up and down for each genotype to get the overall. Then I use set
    terminology to determine which genes are shared and unshared for each
    comparison.

    Args:
        all_degs_pandas (pandas.Data.Frame): Pandas dataframe of
        blueberry genes and DEG comparison columns. Values are Down, Up, or
        No_Change to signify what kind of DEG a gene is.
            Columns:
                Blueberry_Gene: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

        draper_col (str): String of column name to subset out of the
        all_degs_pandas pandas dataframe

        liberty_col (str): String of column name to subset out of the
        all_degs_pandas pandas dataframe

    Returns:
        common_genes (pandas.DataFrame):
                Index:
                    RangeIndex
                Columns:
                    Blueberry_Gene: Object, blueberry gene names
                    str(draper_col): Name varies, Object, values are 'Up' or
                        'Down'
                    str(liberty_col): Name varies, Object, values are 'Up' or
                        'Down'

        unique_draper_genes (pandas.DataFrame):
                Index:
                    RangeIndex
                Columns:
                    Blueberry_Gene: Object, blueberry gene names
                    str(draper_col): Name varies, Object, values are 'Up' or
                        'Down'

        unique_draper_genes (pandas.DataFrame):
                Index:
                    RangeIndex
                Columns:
                    Blueberry_Gene: Object, blueberry gene names
                    str(draper_col): Name varies, Object, values are 'Up' or
                        'Down'


    """
    # Subset the DEG data table by a certain day comparison only for Draper and
    # only keep the genes that are up or down regulated
    all_changes_draper_pandas = all_degs_pandas.loc[
        all_degs_pandas[draper_col] != "No_Change"
    ]
    all_changes_draper_pandas = all_changes_draper_pandas.loc[
        :, [draper_col, "Arabidopsis_Gene"]
    ]
    all_changes_draper_pandas.reset_index(inplace=True)

    # Subset the DEG data table by a certain day comparison only for Liberty and
    # only keep the genes that are up or down regulated
    all_changes_liberty_pandas = all_degs_pandas.loc[
        all_degs_pandas[liberty_col] != "No_Change"
    ]
    all_changes_liberty_pandas = all_changes_liberty_pandas.loc[
        :, [liberty_col, "Arabidopsis_Gene"]
    ]
    all_changes_liberty_pandas.reset_index(inplace=True)

    # The common genes are the INNER merge (i.e the intersection) between the
    # two dataframes
    common_genes = pd.merge(
        all_changes_draper_pandas,
        all_changes_liberty_pandas,
        how="inner",
        on=["Blueberry_Gene", "Arabidopsis_Gene"],
    )

    # Subset the already subsetted Liberty dataset by the genes that are NOT
    # shared with the Draper dataset. This gives a dataset of the uniquely
    # Liberty genes
    unique_liberty_genes = all_changes_liberty_pandas[
        ~all_changes_liberty_pandas["Blueberry_Gene"].isin(
            all_changes_draper_pandas["Blueberry_Gene"]
        )
    ]

    # Subset the already subsetted Draper dataset by the genes that are NOT
    # shared with the Liberty dataset. This gives a dataset of the uniquely
    # Draper genes
    unique_draper_genes = all_changes_draper_pandas[
        ~all_changes_draper_pandas["Blueberry_Gene"].isin(
            all_changes_liberty_pandas["Blueberry_Gene"]
        )
    ]

    # Simple sort of column names
    common_genes.sort_index(axis=1, ascending=True, inplace=True)
    unique_draper_genes = unique_draper_genes.copy(deep=True).sort_index(
        axis=1, ascending=True
    )
    unique_liberty_genes = unique_liberty_genes.copy(deep=True).sort_index(
        axis=1, ascending=True
    )

    # Sort so that the "No_Ortholog" arabidopsis values are at the bottom of
    # the table
    common_genes.sort_values(by=["Arabidopsis_Gene"], inplace=True)
    unique_draper_genes.sort_values(by=["Arabidopsis_Gene"], inplace=True)
    unique_liberty_genes.sort_values(by=["Arabidopsis_Gene"], inplace=True)

    return (
        common_genes,
        unique_draper_genes,
        unique_liberty_genes,
    )


def graph_venn_diagram_deg_per_context(all_degs_pandas, output_dir):
    """
    Create a Venn diagram of the TOTAL DEGs for each genotype at each stage.

    # TODO add the arabidopsis gene column info to the documentation and
    # potentially rename some of these parameters
    Args:
        all_degs_pandas (pandas.Data.Frame): Pandas dataframe of
        blueberry genes and DEG comparison columns. Values are Down, Up, or
        No_Change to signify what kind of DEG a gene is.
            Columns:
                Blueberry_Gene: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

        output_dir (str): Path

    Returns: None, makes graphs and saves them to the output dir. The code also
    saves the pandas DataFrames to disk so that we can have lists of the unique
    and shared genes.
    """
    all_degs_pandas.set_index("Blueberry_Gene", inplace=True)
    columns = sorted(all_degs_pandas.columns.to_list())
    columns.remove("Arabidopsis_Gene")  # NB remove the Arabidopsis gene column
    # to make the below list comprehensions easier
    only_draper_colnames = [genome for genome in columns if "LIB" not in genome]
    only_liberty_colnames = [genome for genome in columns if "DRA" not in genome]

    for draper_col, liberty_col in zip(only_draper_colnames, only_liberty_colnames):
        (
            common_genes,
            unique_draper_genes,
            unique_liberty_genes,
        ) = get_unique_and_non_unique_genes_per_comparison(
            all_degs_pandas.copy(deep=True), draper_col, liberty_col
        )

        count_of_common_genes = len(common_genes)
        count_of_unique_draper_genes = len(unique_draper_genes)
        count_of_unique_liberty_genes = len(unique_liberty_genes)

        day_number = int(draper_col.rstrip("_Direction")[-1])  # MAGIC
        # Magic number just to get the Day integer
        comparison_name = (
            "Day %s: All DEGS of Resistant Genotype (Draper) vs. All DEGs of Susceptible Genotype (Liberty)"
            % day_number
        )

        comparison_filename = "of_Day_" + str(day_number) + ".tsv"

        common_genes.to_csv(
            os.path.join(output_dir, "Shared_DEGS_" + comparison_filename),
            header=True,
            index=False,
            sep="\t",
        )
        unique_draper_genes.to_csv(
            os.path.join(output_dir, "Unique_Draper_DEGS_" + comparison_filename),
            header=True,
            index=False,
            sep="\t",
        )
        unique_liberty_genes.to_csv(
            os.path.join(output_dir, "Unique_Liberty_DEGS_" + comparison_filename),
            header=True,
            index=False,
            sep="\t",
        )

        the_labels = ("Unique Draper DEGs", "Unique Liberty DEGs")
        venn2(
            subsets=(
                count_of_unique_draper_genes,
                count_of_unique_liberty_genes,
                count_of_common_genes,
            ),
            set_colors=("b", "r"),
            set_labels=the_labels,
        )

        # Add another set to make the edges bold. Purely cosmetic
        venn2_circles(
            subsets=(
                count_of_unique_draper_genes,
                count_of_unique_liberty_genes,
                count_of_common_genes,
            )
        )
        plt.title(comparison_name)
        plt.savefig(
            os.path.join(output_dir, "Day" + str(day_number) + "_DEG_VennDiagram.png"),
            bbox_inches="tight",
        )
        plt.clf()


def graph_line_plot_deg_per_context(merged_deg_no_module_ID, output_dir):
    """
    Args:
        merged_deg_no_module_ID (pandas.Data.Frame): Pandas dataframe of
        blueberry genes and DEG comparison columns. Values are Down, Up, or
        No_Change to signify what kind of DEG a gene is.
            Columns:
                Blueberry_Gene: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

        output_dir (str): Path
    """
    all_columns = merged_deg_no_module_ID.columns.to_list()
    all_columns.remove("Blueberry_Gene")
    merged_deg_no_module_ID.set_index("Blueberry_Gene", inplace=True)
    data = merged_deg_no_module_ID.apply(pd.value_counts)
    data.fillna(0, inplace=True)
    columns = sorted(data.columns.to_list())
    only_draper = [genome for genome in columns if "LIB" not in genome]
    only_liberty = [genome for genome in columns if "DRA" not in genome]

    # Draper
    draper_count_dataframe = data[only_draper]
    down_sequence_draper = draper_count_dataframe.loc["Down"].to_list()
    up_sequence_draper = draper_count_dataframe.loc["Up"].to_list()

    # Liberty
    liberty_count_dataframe = data[only_liberty]
    down_sequence_liberty = liberty_count_dataframe.loc["Down"].to_list()
    up_sequence_liberty = liberty_count_dataframe.loc["Up"].to_list()

    x_range = ["C" + str(i) + " vs. " + "T" + str(i) for i in range(1, 8, 1)]
    M = 10
    plt.figure(figsize=(16, 9.5))
    plt.plot(
        x_range,
        up_sequence_draper,
        "^b:",
        label="Up-Regulated in Resistant Genotype (Draper)",
        markersize=M,
    )
    plt.plot(
        x_range,
        down_sequence_draper,
        "vb-",
        label="Down-Regulated in Resistant Genotype (Draper)",
        markersize=M,
    )
    plt.plot(
        x_range,
        up_sequence_liberty,
        "^r:",
        label="Up-Regulated in Susceptible Genotype (Liberty)",
        markersize=M,
    )
    plt.plot(
        x_range,
        down_sequence_liberty,
        "vr-",
        label="Down-Regulated in Susceptible Genotype (Liberty)",
        markersize=M,
    )
    plt.ylabel("Number of Genes")
    plt.xlabel("Comparison Context")
    plt.title("Number of Genes Differentially Expressed Per Context")
    plt.legend()
    plt.savefig(
        os.path.join(output_dir, "Lineplot_number_diffex_genes_per_context.png"),
        bbox_inches="tight",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    # NB input args correspond to Makefile, must be ordered correctly
    parser.add_argument("DEG_dir", type=str, help="output from EdgeR")
    parser.add_argument("GO_ID_Term_file", type=str, help="TODO")
    parser.add_argument(
        "orthology_file",
        type=str,
        help="path to synteny/homology file containing AT BB pairs",
    )
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.DEG_dir = os.path.abspath(args.DEG_dir)
    args.GO_ID_Term_file = os.path.abspath(args.GO_ID_Term_file)
    args.orthology_file = os.path.abspath(args.orthology_file)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE begin reading data
    list_of_deg_paths = read_DEG_dir(args.DEG_dir)
    list_of_deg_pandas = read_all_DEG_direction_files(list_of_deg_paths)
    orthology_table = read_synteny_homology_table(args.orthology_file)

    # NOTE unused while waiting to hear back from Melanie on formatting ideas
    # TODO clean this up if it isn't ultimately used
    # GO_ID_Term_table = read_GO_ID_w_term(args.GO_ID_Term_file)

    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )
    all_degs_and_gene_names = merged_deg_no_module_ID.merge(
        orthology_table, on="Blueberry_Gene", how="left"
    )
    all_degs_and_gene_names.drop(columns=["E_Value", "Point_of_Origin"], inplace=True)
    all_degs_and_gene_names.fillna("No_Ortholog", inplace=True)

    # NOTE this is a unique special result folder
    graph_venn_diagram_deg_per_context(
        all_degs_and_gene_names.copy(deep=True),
        os.path.join(args.output_dir, "Unique_and_Shared_DEGs"),
    )
    graph_line_plot_deg_per_context(
        merged_deg_no_module_ID.copy(deep=True), args.output_dir
    )

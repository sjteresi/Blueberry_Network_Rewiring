#!/usr/bin/env python3

"""
Filter DEG representation data. Determine overlap with modules from
WGCNA.

Create a line plot of DEGs over the time points
Create a Venn Diagram of DEGs over the time points to see overlap between
Liberty and Draper datasets
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from functools import reduce

from src.modules.filter_modules import read_gene_modules_table
from src.module_overlap.module_log2fc_overlap import (
    calc_gene_counts_per_module,
    top_n_modules_and_percent_per_context,
)


def read_DEG_dir(input_directory):
    """
    Reads an input directory and returns a list of file paths that contain
    the magic substring 'Direction.tsv'.

    Args:
        input_directory (str): Path to directory that contains files with our
        magic substring. This directory is parsed to find filenames that match
        the substring.

    Returns:
        deg_files (list): List of absolute file paths
    """
    deg_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("Direction" in f)
    ]  # MAGIC substring for filename recognition

    if len(deg_files) != 28:  # MAGIC number of files, anymore or less will
        # result in error
        raise ValueError(
            """More files were tagged than expected, check to see
                         if files are being tagged appropriately."""
        )

    return deg_files


def read_all_DEG_direction_files(deg_files):
    """
    Reads a list of file paths and reads them as pandas files. Creates a list
    of pandas dataframes. Magic column names are read from the file

    Args:
        deg_files (str): Output from read_log_2fc_dir function. List
        of file paths.

    Returns:
        list_of_read_deg_files (list): List of pandas dataframes
    """
    list_of_read_deg_files = []
    for deg_path in deg_files:
        individual_deg_file = pd.read_csv(deg_path, sep="\t", header="infer")

        if len(os.path.splitext(os.path.basename(deg_path))) > 2:
            # MAGIC remove filename extension, does not work if multiple periods
            # are in filename.
            raise ValueError(
                """There were multiple periods in your file path,
                             this will cause an error."""
            )
        analysis_context = os.path.splitext(os.path.basename(deg_path))[0]
        analysis_context = analysis_context.rstrip("_Direction")  # NB strip
        # the string from the end so the groupby unstack command in
        # convert_direction_integer_to_string gives us a more appropriate
        # column name

        individual_deg_file.rename(
            columns={
                "Direction_Differentially_Regulated": analysis_context,
                "Gene_Name": "Blueberry_Gene",
            },
            inplace=True,
        )  # MAGIC column name
        individual_deg_file = convert_direction_integer_to_string(
            individual_deg_file, analysis_context
        )
        list_of_read_deg_files.append(individual_deg_file)

    return list_of_read_deg_files


def panda_info(panda_file):
    print(panda_file)
    print(panda_file.shape)
    print(panda_file.columns)
    print(panda_file.info())


def convert_direction_integer_to_string(deg_panda_file, analysis_context):
    """
    Take the output from EdgeR (the direction files), where a gene can either
    be 0, 1, or -1 and convert those integer files to strings signifying the
    direction of regulation (i.e up, down, or no change)

    Args:
        deg_panda_file (pandas.DataFrame):
            Index:
                RangeIndex
            Columns:
                Blueberry_Gene: Object, dtype: object
                analysis_context: Integer, dtype: int64
            Shape:
                (93583, 2)

        analysis_context (string):
            String representing the column name and comparison context
            of the EdgeR data.

    Returns: deg_panda_file (pandas.DataFrame):
        Index:
            RangeIndex
        Columns:
            Blueberry_Gene: Object, dtype: object
            analysis_context + _'Direction': Object, dtype: object
        Shape:
            (93583, 2)
    """
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == -1, analysis_context + "_Direction"
    ] = "Down"
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == 1, analysis_context + "_Direction"
    ] = "Up"
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == 0, analysis_context + "_Direction"
    ] = "No_Change"
    deg_panda_file.drop(columns=analysis_context, inplace=True)

    return deg_panda_file


def count_direction_of_regulation_per_module_deg(all_merged):
    """
    Converts the Log2FC columns of all_merged into counts per module. Index is
    now the module ID not the gene. There is a count of genes in each context

    Args:
        all_merged (pandas.DataFrame):
            Index:
                RangeIndex:
            Columns:
                Blueberry_Gene: Object, dtype: object
                Number_of_Blueberry_Genes_in_Module: Integer, dtype: int64
                Module_Color: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

    Returns:
        counts (pandas.DataFrame):
            Index:
                Module_Color: Length 89: Object, dtype: object

            Columns:
                Expression Contexts (MAGIC VARIES): Object: dtype: object
    """
    unwanted_columns = [
        "Blueberry_Gene",
        "Module_Color",
        "Number_of_Blueberry_Genes_in_Module",
    ]
    columns_of_direction = all_merged.columns.to_list()
    columns_of_direction = [
        column_name
        for column_name in columns_of_direction
        if column_name not in unwanted_columns
    ]

    counts_per_module_comparison = []
    for direction_column in columns_of_direction:
        x = (
            all_merged.groupby(by="Module_Color", dropna=True)[direction_column]
            .value_counts()
            .unstack()
        )
        comparison_name = x.columns.name
        x = x.add_prefix(comparison_name + "_")
        x = x.rename_axis(None, axis=1)
        counts_per_module_comparison.append(x)
        col_list = [col for col in x.columns if "No_Change" not in col]
        x[direction_column + "_Total_DEGS"] = x[col_list].sum(axis=1)

    # Perform a pandas merge on a list
    counts = reduce(
        lambda left, right: pd.merge(left, right, on="Module_Color", how="outer"),
        counts_per_module_comparison,
    )
    return counts


def get_counts_of_unique_and_non_unique_genes_per_comparison(
    merged_deg_no_module_ID, draper_col, liberty_col
):
    """
    Gets the counts of unique and non-unique DEGs per genotype comparison.
    Basically does Draper C1 and T1 vs Liberty C1 and T1. I am grouping the
    DEGs of up and down for each genotype to get the overall. Then I use set
    terminology to determine which genes are shared and unshared for each
    comparison.

    Args:
        merged_deg_no_module_ID (pandas.Data.Frame): Pandas dataframe of
        blueberry genes and DEG comparison columns. Values are Down, Up, or
        No_Change to signify what kind of DEG a gene is.
            Columns:
                Blueberry_Gene: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

        draper_col (str): String of column name to subset out of the
        merged_deg_no_module_ID pandas dataframe

        liberty_col (str): String of column name to subset out of the
        merged_deg_no_module_ID pandas dataframe

    Returns:
        awful tuple lol please refactor if using in the future:
            count_of_common_genes (int):
            count_of_unique_draper_genes (int)
            count_of_unique_liberty_genes (int)
    """
    # Get smaller panda dataframes with only the relevant data
    all_changes_draper_series = merged_deg_no_module_ID[draper_col].loc[
        merged_deg_no_module_ID[draper_col] != "No_Change"
    ]
    all_changes_liberty_series = merged_deg_no_module_ID[liberty_col].loc[
        merged_deg_no_module_ID[liberty_col] != "No_Change"
    ]
    common_genes = set(all_changes_draper_series.index).intersection(
        set(all_changes_liberty_series.index)
    )
    count_of_common_genes = len(common_genes)
    count_of_unique_draper_genes = (
        len(all_changes_draper_series) - count_of_common_genes
    )
    count_of_unique_liberty_genes = (
        len(all_changes_liberty_series) - count_of_common_genes
    )
    return (
        count_of_common_genes,
        count_of_unique_draper_genes,
        count_of_unique_liberty_genes,
    )


def graph_venn_diagram_deg_per_context(merged_deg_no_module_ID, output_dir):
    """
    Create a Venn diagram of the TOTAL DEGs for each genotype at each stage.

    Args:
        merged_deg_no_module_ID (pandas.Data.Frame): Pandas dataframe of
        blueberry genes and DEG comparison columns. Values are Down, Up, or
        No_Change to signify what kind of DEG a gene is.
            Columns:
                Blueberry_Gene: Object, dtype: object
                Expression Contexts (MAGIC VARIES): Object: dtype: object

        output_dir (str): Path

    Returns: None, just makes graphs and saves them to the output dir.
    """
    all_columns = merged_deg_no_module_ID.columns.to_list()
    all_columns.remove("Blueberry_Gene")
    merged_deg_no_module_ID.set_index("Blueberry_Gene", inplace=True)
    columns = sorted(merged_deg_no_module_ID.columns.to_list())
    only_draper_colnames = [genome for genome in columns if "LIB" not in genome]
    only_liberty_colnames = [genome for genome in columns if "DRA" not in genome]

    for draper_col, liberty_col in zip(only_draper_colnames, only_liberty_colnames):
        (
            count_of_common_genes,
            count_of_unique_draper_genes,
            count_of_unique_liberty_genes,
        ) = get_counts_of_unique_and_non_unique_genes_per_comparison(
            merged_deg_no_module_ID, draper_col, liberty_col
        )
        day_number = int(draper_col.rstrip("_Direction")[-1])  # MAGIC
        # Magic number just to get the Day integer
        comparison_name = (
            "Day %s: All DEGS of Resistant Genotype (Draper) vs. All DEGs of Susceptible Genotype (Liberty)"
            % day_number
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
        os.path.join(output_dir, "lineplot_number_diffex_genes_per_context.png"),
        bbox_inches="tight",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("DEG_dir", type=str, help="output from EdgeR")
    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.DEG_dir = os.path.abspath(args.DEG_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE read the input datasets
    blueberry_genes_and_module_colors_table = read_gene_modules_table(
        args.genes_and_module_colors
    )
    gene_counts_per_module = calc_gene_counts_per_module(
        blueberry_genes_and_module_colors_table
    )

    # NOTE begin the analysis
    list_of_deg_paths = read_DEG_dir(args.DEG_dir)
    list_of_deg_pandas = read_all_DEG_direction_files(list_of_deg_paths)

    merged_deg_no_module_ID = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        list_of_deg_pandas,
    )

    graph_venn_diagram_deg_per_context(merged_deg_no_module_ID, args.output_dir)
    graph_line_plot_deg_per_context(merged_deg_no_module_ID, args.output_dir)

    merged_deg_w_module_ID = pd.merge(
        merged_deg_no_module_ID,
        blueberry_genes_and_module_colors_table,
        on="Blueberry_Gene",
        how="outer",
    )
    degs_and_modules = merged_deg_w_module_ID.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )

    modules_and_count_per_context = count_direction_of_regulation_per_module_deg(
        degs_and_modules,
    )
    # Drop the 'No_Change' columns
    list_of_no_change_cols = [
        f for f in modules_and_count_per_context.columns if "No_Change" in f
    ]
    modules_and_count_per_context.drop(columns=list_of_no_change_cols, inplace=True)

    # Get the gene counts per module information again because it was lost in
    # the previous function.
    modules_and_count_per_context = modules_and_count_per_context.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )
    modules_and_count_per_context.set_index("Module_Color", inplace=True)

    # Divide to get the percentage
    percentages = modules_and_count_per_context.apply(
        lambda row: row / row["Number_of_Blueberry_Genes_in_Module"], axis=1
    )  # NB calculate percent, make it a new pandas dataframe
    percentages.drop(
        columns=["Number_of_Blueberry_Genes_in_Module"], inplace=True
    )  # NB remove column
    percentages = percentages.add_prefix("Percent_")  # NB add prefix to column

    # NOTE FUTURE
    # Get gene counts per module again because the number of blueberry genes
    # per module column was set to 1 because of division, in future figure out
    # how to do if else statements in pandas apply lambda. Not sure how to keep
    # it from acting on a specific column, would I first need to provide a
    # mask instead of doing if else?
    percentages = percentages.merge(
        gene_counts_per_module, on="Module_Color", how="outer"
    )
    percentages.fillna(np.nan, inplace=True)
    percentages.set_index("Module_Color", inplace=True)
    top_n_modules_and_percent_per_context(percentages, args.output_dir)

    # MAGIC filename
    percentages.to_csv(
        os.path.join(args.output_dir, "DEG_percentages_in_modules.tsv"),
        sep="\t",
        header=True,
        index=True,
    )

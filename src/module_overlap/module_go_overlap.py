#!/usr/bin/env python3

"""
Filter GO data. Determine overlap with modules from WGCNA.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os
from functools import reduce

from src.read_tables_and_dir import (
    read_gene_modules_table,
    read_interesting_GO_terms,
    read_FPKM_or_TPM,
    read_synteny_homology_table,
    read_GO_ID_w_term,
    read_GO_dir,
)
from src.module_overlap.module_log2fc_overlap import (
    calc_gene_counts_per_module,
    read_log_2fc_dir,
    read_all_log_2fc_files,
)


def read_GO_files_and_merge(go_files):
    """
    Reads multiple GO files (output from TOPGO saying if a GO term was found to
    be enriched in a given module), and concatenates and merge them into a useful
    format for downstream analysis.

    Return a list of pandas dataframes, one pandaframe for each module ID. Each
    pandaframe is the GO terms and IDs that were recovered in that module

    Args:
        go_files (str): Output from read_log_2fc_dir function. List
        of file paths.

    Returns:
        list_pandas_per_module_w_go (list of pandas.DataFrame): Each pandas dataframe has form:
            Index:
                RangeIndex:
            Columns:
                GO_ID: Object, dtype: object
                GO_TERM: Object, dtype: object
                GO_Term_Type: Object, dtype: object
                Module_Name (VARIES): Object, dtype: object
    """
    dictionary_of_module_w_panda_list = {}
    for go_file_path in go_files:
        individual_go_panda = pd.read_csv(go_file_path, sep="\t", header="infer")
        individual_go_panda.rename(
            columns={
                "GO.ID": "GO_ID",
                "Term": "GO_TERM",
            },
            inplace=True,
        )  # MAGIC column name
        individual_go_panda.drop(
            columns=[
                "Annotated",
                "Significant",
                "Expected",
                "classicFisher",
                "Overrepresented",
            ],
            inplace=True,
        )

        if len(os.path.splitext(os.path.basename(go_file_path))) > 2:
            # MAGIC remove filename extension, does not work if multiple periods
            # are in filename.
            raise ValueError(
                """There were multiple periods in your file path,
                             this will cause an error."""
            )
        # MAGIC to get filename
        analysis_context = os.path.splitext(os.path.basename(go_file_path))[0]

        # MAGIC string split on underscores in filename, there is a common
        # structure to the filename that allows me to do this
        module_name = analysis_context.split("_")[0]
        analysis_type = analysis_context.split("_")[4]

        individual_go_panda["GO_Term_Type"] = analysis_type
        individual_go_panda[module_name] = "Recovered"

        if module_name not in dictionary_of_module_w_panda_list:
            dictionary_of_module_w_panda_list[module_name] = []
        dictionary_of_module_w_panda_list[module_name].append(individual_go_panda)

    # NB combine the dataframes
    list_pandas_per_module_w_go = [
        pd.concat(val) for key, val in dictionary_of_module_w_panda_list.items()
    ]

    # Perform a pandas merge on a list
    all_go_merged = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=["GO_ID", "GO_TERM", "GO_Term_Type"],
            how="outer",
        ),
        list_pandas_per_module_w_go,
    )
    return all_go_merged


def gen_melanie_table(interesting_go_terms):
    """
    Takes the incoming panda frame, which has GO terms as rows and modules as
    columns, and re-orients to have modules as rows ad go terms in the columns.
    It shows the GO terms and IDs of interest as they appear in the modules.
    This is arguably a little easier to read, there are GO terms that repeat but
    that is on purpose.

    Args:
        interesting_go_terms (pandas.DataFrame):

    Returns:
        data_table (pandas.DataFrame):
    """

    # The list comprehension creates a whitelist for us to easily perform
    # calculations on the relevant columns
    # NOTE errors will happen in future work if column names are not added to
    # the below list comprehension.
    module_names = [
        module_name
        for module_name in interesting_go_terms.columns.to_list()
        if module_name not in ["GO_Term_Type", "GO_TERM", "GO_ID"]
    ]

    to_concat = []
    for module_name in module_names:
        mini_df = interesting_go_terms.loc[
            interesting_go_terms[module_name] == "Recovered"
        ]
        mini_df = mini_df.loc[:, ["GO_Term_Type", "GO_TERM", "GO_ID", module_name]]
        col = mini_df.pop(module_name)
        mini_df.insert(0, col.name, col)  # MAGIC, put in first column
        mini_df.rename(columns={module_name: "Module_ID"}, inplace=True)
        mini_df["Module_ID"] = col.name
        to_concat.append(mini_df)

    data_table = pd.concat(to_concat)
    data_table.sort_values(by=["GO_ID"], inplace=True)
    return data_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("go_dir", type=str, help="output from topGO")
    parser.add_argument(
        "genes_and_module_colors",
        type=str,
        help="file representing the blueberry genes and their module IDs",
    )
    parser.add_argument(
        "interesting_go_terms",
        type=str,
        help="file representing the interesting GO terms received from Melanie",
    )
    parser.add_argument(
        "log2fc_expression_data",
        type=str,
        help="file representing the gene expression data",
    )
    parser.add_argument(
        "orthologs",
        type=str,
        help="file representing the orthlog data",
    )
    parser.add_argument(
        "GO_ID_Term_file",
        type=str,
        help="""filtered GO file
                        from custom script""",
    )
    parser.add_argument("output_dir", type=str, help="output_dir path")

    # NB set args
    args = parser.parse_args()
    args.go_dir = os.path.abspath(args.go_dir)
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.interesting_go_terms = os.path.abspath(args.interesting_go_terms)
    args.log2fc_expression_data = os.path.abspath(args.log2fc_expression_data)
    args.orthologs = os.path.abspath(args.orthologs)
    args.GO_ID_Term_file = os.path.abspath(args.GO_ID_Term_file)
    args.output_dir = os.path.abspath(args.output_dir)

    # NB set logging
    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NB
    go_file_paths = read_GO_dir(args.go_dir)
    all_go_merged = read_GO_files_and_merge(go_file_paths)

    # NOTE
    # all_go_merged is the entire raw data. We must now begin to filter it
    # Remove rows where the entire row is NA
    # Remove columns (modules) where the entire module is NA
    all_go_merged.dropna(axis=0, how="all", inplace=True)
    all_go_merged.dropna(axis=1, how="all", inplace=True)
    all_go_merged.fillna("Not_Recovered", inplace=True)
    # NB save for future even though it may not be all that usful
    all_go_merged.to_csv(
        os.path.join(args.output_dir, "TopGO_All_Merged_NoFilter.tsv"),
        sep="\t",
        header=True,
    )

    # Melanie built a dataset of terms she already had an eye on.
    # We will use this to subset
    interesting_go_terms = read_interesting_GO_terms(args.interesting_go_terms)
    all_go_merge_interesting_subset = all_go_merged[
        all_go_merged["GO_ID"].isin(interesting_go_terms["GO_ID"])
    ]
    # NB save for future, this is the all_go_merged data subset by the
    # interesting genes
    all_go_merge_interesting_subset.to_csv(
        os.path.join(args.output_dir, "Interesting_GO_Term_Module_Representation.tsv"),
        sep="\t",
        header=True,
        index=False,
    )

    # NOTE
    # Read additional datasets
    # Expression
    list_log_2fc_files = read_log_2fc_dir(args.log2fc_expression_data)
    pandas_log_2fc_files = read_all_log_2fc_files(list_log_2fc_files)

    # NOTE begin the analysis
    # NB do a pandas merge over a list, merge all the log2fc files (which have
    # different columns), but a common blueberry gene column
    expression_data = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        pandas_log_2fc_files,
    )
    expression_data.fillna("NA", inplace=True)

    ortholog_data = read_synteny_homology_table(args.orthologs)
    GO_ID_Term_table = read_GO_ID_w_term(args.GO_ID_Term_file)
    GO_ID_Term_table.drop(columns=["GO_Term_Description"], inplace=True)
    # Now begin to add additional information to the table
    orthos_and_expression = pd.merge(
        expression_data, ortholog_data, how="outer", on="Blueberry_Gene"
    )

    data = gen_melanie_table(all_go_merge_interesting_subset)

    # Add the arabidopsis genes to the table
    data = data.merge(GO_ID_Term_table, how="inner", on="GO_ID")

    # Add the blueberry genes and gene expression columns to the table
    data = data.merge(orthos_and_expression, how="inner", on="Arabidopsis_Gene")
    data.drop(columns=["E_Value", "Point_of_Origin", "GO_Term_Type"], inplace=True)

    # Save the data table
    data.to_csv(
        os.path.join(
            args.output_dir, "Interesting_GO_Terms_Module_Representation_COMPLETE.tsv"
        ),
        sep="\t",
        header=True,
        index=False,
    )

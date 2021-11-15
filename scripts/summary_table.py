#!/usr/bin/env python3

"""
Create the union of dataframes of differentially expressed genes, Arabidopsis
ortholog, Arabidopsis GO term, and blueberry gene network module identity.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os
from functools import reduce

from scripts.modules.filter_modules import read_synteny_homology_table
from scripts.modules.filter_modules import read_gene_modules_table


def read_TPM_expression_table(input_file):
    data = pd.read_csv(input_file, sep="\t", header="infer")
    data.rename(columns={"Gene_Name": "Blueberry_Gene"}, inplace=True)
    return data


def read_diff_ex_direction_directory(input_directory):
    only_direction_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("_Summary" not in f)
    ]
    comparison_groupings = [
        os.path.basename(os.path.splitext(f)[0]).rstrip("_Direction")
        for f in only_direction_files
    ]

    # NB sort alphabetically so we can zip and then pair up
    only_direction_files.sort()
    comparison_groupings.sort()

    table_list = []
    for direction_data, comparison_group in zip(
        only_direction_files, comparison_groupings
    ):
        table = read_diff_ex_direction_file(direction_data, comparison_group)
        table_list.append(table)
    return table_list


def read_diff_ex_direction_file(input_file, comparison_group):
    data = pd.read_csv(input_file, sep="\t", header="infer")
    data.rename(
        columns={
            "Gene_Name": "Blueberry_Gene",
            "Direction_Differentially_Regulated": comparison_group,
        },
        inplace=True,
    )
    return data


def read_arabidopsis_go_table(input_file):
    data = pd.read_csv(
        input_file, sep="\t", header=None, names=["Arabidopsis_Gene", "GO_Terms"]
    )
    data["GO_Terms"] = [x.replace(" ", "").split(",") for x in data["GO_Terms"]]
    return data


# TODO consider adding TPM/FPKM?

# TODO create function to merge all the datasets
def merge_dataframes(
    synteny_homology_data,
    arabidopsis_go_data,
    blueberry_genes_and_module_colors,
    blueberry_tpm,
    diff_ex_panda_list,
):

    # NOTE merge blueberry gene / module color with blueberry gene / diff ex
    # status
    all_diff_ex_merged = reduce(
        lambda left, right: pd.merge(left, right, on="Blueberry_Gene", how="outer"),
        diff_ex_panda_list,
    )
    merged = pd.merge(
        all_diff_ex_merged,
        blueberry_genes_and_module_colors,
        on="Blueberry_Gene",
        how="outer",
    )

    # NOTE merge in synteny homology data to get Arabidopsis gene, some
    # Arabidopsis genes will repeat, and there are some blueberry genes from
    # scaffolds that aren't the main chromosomes. They will have NAs in the
    # diff ex and module set.
    # This has all of the blueberry genes, no duplicates
    # SQL full outer join
    merged = pd.merge(merged, synteny_homology_data, on="Blueberry_Gene", how="outer")

    # NOTE merge in blueberry TPM data
    merged = pd.merge(merged, blueberry_tpm, on="Blueberry_Gene", how="outer")

    # NOTE get the GO data for an Arabidopsis gene only if it is already in the
    # data table, SQL left outer join
    complete = pd.merge(merged, arabidopsis_go_data, on="Arabidopsis_Gene", how="left")

    return complete


# TODO create function to save final output or perform one-liner in main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # TODO give description str
    parser.add_argument("genes_and_module_colors", type=str, help="")

    # TODO give description str
    parser.add_argument("synteny_homology_table", type=str, help="")

    # NOTE 'direction' files will have to be iterated over
    # TODO give description str
    parser.add_argument("diff_ex_dir", type=str, help="")

    # TODO give description str
    parser.add_argument("arabidopsis_go_terms", type=str, help="")

    # TODO give description str
    parser.add_argument("blueberry_tpm", type=str, help="")

    parser.add_argument("output_dir", type=str, help="parent path to output results")

    args = parser.parse_args()
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.synteny_homology_table = os.path.abspath(args.synteny_homology_table)
    args.diff_ex_dir = os.path.abspath(args.diff_ex_dir)
    args.arabidopsis_go_terms = os.path.abspath(args.arabidopsis_go_terms)
    args.blueberry_tpm = os.path.abspath(args.blueberry_tpm)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    print()

    synteny_homology_data = read_synteny_homology_table(args.synteny_homology_table)
    arabidopsis_go_data = read_arabidopsis_go_table(args.arabidopsis_go_terms)
    blueberry_genes_and_module_colors = read_gene_modules_table(
        args.genes_and_module_colors
    )
    blueberry_tpm = read_TPM_expression_table(args.blueberry_tpm)
    diff_ex_panda_list = read_diff_ex_direction_directory(args.diff_ex_dir)
    # NB all primary datafiles now read in

    complete = merge_dataframes(
        synteny_homology_data,
        arabidopsis_go_data,
        blueberry_genes_and_module_colors,
        blueberry_tpm,
        diff_ex_panda_list,
    )
    # TODO get filename to save
    complete.to_csv("test.tsv", sep="\t", header=True, index=False)

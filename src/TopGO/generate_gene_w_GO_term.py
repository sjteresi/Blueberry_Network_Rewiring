#!/usr/bin/env python3

"""
Filter a GO term database of Arabidopsis genes.
Generate two dataframes:
    Arabidopsis genes as rows, GO IDs as weird list format that TopGO wants

    GO IDs as rows, GO terms (string) in second column
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os


def generate_go_id_and_term_table(go_master_file):
    """
    TODO
    """
    column_names = [
        "locus_name",
        "TAIR_accession",
        "object_name",
        "relationship_type",
        "GO_term",
        "GO_ID",
        "TAIR_keyword_ID",
        "aspect",
        "GOslim_term",
        "evidence_code",
        "evidence_description",
        "evidence_with",
        "reference",
        "annotator",
        "date",
    ]
    columns_to_use = ["GO_ID", "GO_term"]
    data = pd.read_csv(
        go_master_file,
        sep="\t",
        comment="!",
        header=None,
        names=column_names,
        usecols=columns_to_use,
    )
    data.drop_duplicates(inplace=True)
    data.rename(columns={"GO_term": "GO_Term"}, inplace=True)
    return data


def read_master_GO_table(go_master_file):
    """
    Take a master GO table from TAIR and preprocess to remove gene IDs we do
    not want.

    Args:
        go_master_file (str to filepath):
    """
    column_names = [
        "locus_name",
        "TAIR_accession",
        "object_name",
        "relationship_type",
        "GO_term",
        "GO_ID",
        "TAIR_keyword_ID",
        "aspect",
        "GOslim_term",
        "evidence_code",
        "evidence_description",
        "evidence_with",
        "reference",
        "annotator",
        "date",
    ]
    columns_to_use = ["locus_name", "GO_ID"]
    data = pd.read_csv(
        go_master_file,
        sep="\t",
        comment="!",
        header=None,
        names=column_names,
        usecols=columns_to_use,
    )

    # Now we need to group by locus_name (rename to gene_name), and condense
    # unique GO_IDs to one gene row.
    data.rename(columns={"locus_name": "gene_name"}, inplace=True)

    # use the group by to gather all GO terms for one gene ID
    data = data.groupby("gene_name")["GO_ID"].apply(",".join).reset_index()

    # Convert the big string of GO terms into a list
    data["GO_ID"] = data["GO_ID"].apply(lambda x: x.split(","))

    # convert list to set to remove duplicate GO terms
    data["GO_ID"] = data["GO_ID"].apply(lambda x: set(x))

    # get in format for TOPGO: Gene <TAB> list of GO terms
    data["GO_ID"] = data["GO_ID"].apply(", ".join)

    # Remove any non-'AT' genes.
    # MAGIC, regex to get AT (1-5) G genes. Don't want mitochondrial or
    # chloroplast or uniquely named genes.
    data = data[data["gene_name"].str.contains("^AT[1-5]G", regex=True)]
    data = data.sort_values(by=["gene_name"])

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("go_master_file", type=str, help="downloaded from TAIR")
    parser.add_argument(
        "go_output_path", type=str, help="output path of filtered GO file"
    )
    args = parser.parse_args()

    args.go_master_file = os.path.abspath(args.go_master_file)
    args.go_output_path = os.path.abspath(args.go_output_path)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    data = read_master_GO_table(args.go_master_file)
    data.to_csv(
        os.path.join(args.go_output_path, "ArabidopsisGene_w_GO.tsv"),
        sep="\t",
        index=False,
        header=False,
    )

    go_id_table = generate_go_id_and_term_table(args.go_master_file)
    go_id_table.to_csv(
        os.path.join(args.go_output_path, "GO_ID_w_Term.tsv"),
        sep="\t",
        index=False,
        header=False,
    )

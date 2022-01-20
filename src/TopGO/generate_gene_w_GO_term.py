#!/usr/bin/env python3

"""
Filter a GO term database of Arabidopsis genes.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os


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
    parser.add_argument("output_dir", type=str, help="parent path to output results")
    args = parser.parse_args()

    args.go_master_file = os.path.abspath(args.go_master_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    data = read_master_GO_table(args.go_master_file)

    # MAGIC filename
    file_output_name = os.path.abspath(
        os.path.join(args.output_dir, "ArabidopsisGene_w_GO.tsv")
    )
    data.to_csv(file_output_name, sep="\t", index=False, header=False)

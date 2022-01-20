#!/usr/bin/env python3

"""
Generate an FPKM table of a blueberry gene and its Arabidopsis syntelog
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import os


def read_fpkm_table(fpkm_table_path):
    fpkm_table = pd.read_csv(fpkm_table_path, sep="\t", header="infer")
    fpkm_table.rename(columns={"Gene_Name": "Blueberry_Gene"}, inplace=True)
    fpkm_table.set_index("Blueberry_Gene", inplace=True)
    return fpkm_table


def read_syntelog_table(syntelog_table_path):
    syntelog_table = pd.read_csv(
        syntelog_table_path, sep="\t", header="infer", index_col="Blueberry_Gene"
    )
    return syntelog_table


def combine_tables(fpkm_table, syntelog_table):
    combined_table = fpkm_table.join(syntelog_table)
    combined_table.fillna("Not_Available", inplace=True)
    return combined_table


def save_output_table(output_table, output_path, output_filename_string, logger):
    output_file_path = os.path.join(output_path, output_filename_string)
    logger.info("Saving file to: %s" % output_file_path)
    output_table.to_csv(output_file_path, sep="\t", header=True, index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("fpkm_table_path", type=str, help="parent path of fpkm table")
    parser.add_argument(
        "syntelog_table_path", type=str, help="parent path of syntelog table"
    )
    parser.add_argument("output_dir", type=str, help="parent path to output results")
    args = parser.parse_args()
    args.fpkm_table_path = os.path.abspath(args.fpkm_table_path)
    args.syntelog_table_path = os.path.abspath(args.syntelog_table_path)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    fpkm_table = read_fpkm_table(args.fpkm_table_path)
    syntelog_table = read_syntelog_table(args.syntelog_table_path)
    combined_table = combine_tables(fpkm_table, syntelog_table)
    save_output_table(
        combined_table, args.output_dir, "Expression_Table_w_Syntelogs.tsv", logger
    )

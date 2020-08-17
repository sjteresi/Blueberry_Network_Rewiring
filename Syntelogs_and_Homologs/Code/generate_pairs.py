#!/usr/bin/env python3

"""
Master code file. Control filtration of syntelog data and homolog data. Manage
merging of the two dataframes and orchestrate all file commands.
"""

__author__ = "Scott Teresi"

import argparse
import os

import logging
import coloredlogs
import pandas as pd

from Code.import_syntelogs import import_syntelogs
from Code.import_homologs import import_homologs
from Code.syntelog_data import Syntelog_Data
from Code.homolog_data import Homolog_Data
from Code.merge_homo_synt import merge_homo_synt
from Code.merged_all_data import Merged_Data


def process(syntelog_input_file, homolog_input_file, data_output_path):
    # Import the data from raw file
    logger.info("Working on syntelogs...")
    syntelogs = import_syntelogs(syntelog_input_file)

    # Wrap the data
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    instance_Syntelog_Data.save_to_disk(
        os.path.join(data_output_path, "set_syntelogs.tsv")
    )

    # Import the data from raw file
    logger.info("Working on homologs...")
    homologs = import_homologs(homolog_input_file)

    # Wrap the data
    instance_Homolog_Data = Homolog_Data(homologs)
    # Save to disk
    instance_Homolog_Data.save_to_disk(
        os.path.join(data_output_path, "set_homologs.tsv")
    )

    # Merge the data
    logger.info("Merging the data...")
    merged_all = merge_homo_synt(instance_Syntelog_Data, instance_Homolog_Data)

    # Wrap the data
    instance_Merged_Data = Merged_Data(merged_all)
    # Save to disk
    instance_Merged_Data.save_to_disk(
        os.path.join(data_output_path, "merged_homo_and_syn.tsv")
    )


if __name__ == "__main__":
    """Command line interface to link syntelogs together."""

    parser = argparse.ArgumentParser(description="Filter syntelogs")
    path_main = os.path.abspath(__file__)
    parser.add_argument(
        "syntelog_input_file", type=str, help="parent path of syntelog file"
    )
    parser.add_argument(
        "homolog_input_file", type=str, help="parent path of homolog file"
    )
    parser.add_argument(
        "--output_directory",
        type=str,
        help="parent path of output directory",
        default=os.path.join(path_main, "../../data_output"),
    )

    parser.add_argument(
        "--input_directory",
        type=str,
        help="parent path of input directory",
        default=os.path.join(path_main, "../../data_input"),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.homolog_input_file = os.path.abspath(args.homolog_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    args.input_directory = os.path.abspath(args.input_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    logger.info("Starting filtration...")
    process(args.syntelog_input_file, args.homolog_input_file, args.output_directory)

#!/usr/bin/env python3

"""
Filter syntelogs from CoGe output file
"""

__author__ = "Scott Teresi"

import argparse
import os

import logging
import coloredlogs
import pandas as pd

from Syntelogs.import_syntelogs import import_syntelogs
from Syntelogs.syntelog_data import Syntelog_Data


def process(input_gene_data, data_output_path):
    # Import the data from raw file
    syntelogs = import_syntelogs(input_gene_data)

    # Wrap the data
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    instance_Syntelog_Data.save_shortened_to_disk(
        os.path.join(data_output_path, "save_test.tsv")
    )


if __name__ == "__main__":
    """Command line interface to link syntelogs together."""

    parser = argparse.ArgumentParser(description="Filter syntelogs")
    path_main = os.path.abspath(__file__)
    parser.add_argument("genes_input_file", type=str, help="parent path of gene file")
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
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    args.input_directory = os.path.abspath(args.input_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    process(args.genes_input_file, args.output_directory)

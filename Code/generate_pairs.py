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

from Code.import_syntelogs import import_syntelogs
from Code.import_orthologs import import_orthologs
from Code.syntelog_data import Syntelog_Data
from Code.ortholog_data import Ortholog_Data
from Code.merge_ortho_synt import merge_ortho_synt
from Code.merged_all_data import Merged_Data


def process(syntelog_input_file, ortholog_input_file, data_output_path):
    # Import the data from raw file
    syntelogs = import_syntelogs(syntelog_input_file)

    # Wrap the data
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    instance_Syntelog_Data.save_to_disk(
        os.path.join(data_output_path, "set_syntelogs.tsv")
    )

    # Import the data from raw file
    orthologs = import_orthologs(ortholog_input_file)

    # Wrap the data
    instance_Ortholog_Data = Ortholog_Data(orthologs)
    instance_Ortholog_Data.save_to_disk(
        os.path.join(data_output_path, "set_orthologs.tsv")
    )

    # Merge the data
    merged_all = merge_ortho_synt(instance_Syntelog_Data, instance_Ortholog_Data)

    # Wrap the data
    instance_Merged_Data = Merged_Data(merged_all)
    instance_Merged_Data.save_to_disk(
        os.path.join(data_output_path, "merged_orth_and_syn.tsv")
    )


if __name__ == "__main__":
    """Command line interface to link syntelogs together."""

    parser = argparse.ArgumentParser(description="Filter syntelogs")
    path_main = os.path.abspath(__file__)
    parser.add_argument(
        "syntelog_input_file", type=str, help="parent path of syntelog file"
    )
    parser.add_argument(
        "ortholog_input_file", type=str, help="parent path of ortholog file"
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
    args.ortholog_input_file = os.path.abspath(args.ortholog_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    args.input_directory = os.path.abspath(args.input_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    process(args.syntelog_input_file, args.ortholog_input_file, args.output_directory)

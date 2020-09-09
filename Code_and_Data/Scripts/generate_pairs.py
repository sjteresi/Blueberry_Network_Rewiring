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

from Code_and_Data.Scripts.import_syntelogs import import_syntelogs
from Code_and_Data.Scripts.import_syntelogs import Syntelog_Data
from Code_and_Data.Scripts.import_homologs import import_homologs
from Code_and_Data.Scripts.import_homologs import Homolog_Data
from Code_and_Data.Scripts.merge_homo_synt import merge_homo_synt
from Code_and_Data.Scripts.merge_homo_synt import Merged_Data
from Code_and_Data.Scripts.import_diff_exp import ExpData
from Code_and_Data.Scripts.union_data import Union_Data


def process(
    syntelog_input_file, homolog_input_file, data_output_path, diff_exp_dir, stat_type,
):
    # Import the synteny data from raw file
    logger.info("Working on syntelogs...")
    syntelogs = import_syntelogs(syntelog_input_file)

    # Wrap the data
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    instance_Syntelog_Data.save_to_disk(
        os.path.join(data_output_path, "set_syntelogs.tsv")
    )

    # Import the homology (BLAST) data from raw file
    logger.info("Working on homologs...")
    homologs = import_homologs(homolog_input_file)

    # Wrap the data
    instance_Homolog_Data = Homolog_Data(homologs)
    # Save to disk
    instance_Homolog_Data.save_to_disk(
        os.path.join(data_output_path, "set_homologs.tsv")
    )

    # Merge the synteny and homology data
    logger.info("Merging the data...")
    merged_all = merge_homo_synt(instance_Syntelog_Data, instance_Homolog_Data)

    # Wrap the data
    instance_Merged_Data = Merged_Data(merged_all)
    # Save to disk
    instance_Merged_Data.save_to_disk(
        os.path.join(data_output_path, "merged_homo_and_syn.tsv")
    )
    # synteny_count, blast_count = instance_Merged_Data.total_count_summary

    # Get differential expression data
    logger.info("Getting the differential expression data...")
    instance_Exp_Data = ExpData(diff_exp_dir)

    # Get union of blueberry genes between diff exp and synteny/homology,
    # and thus the corresponding Arabidopsis gene
    logger.info("Find the union between the diff exp data and the merged data")
    union_data_output_dir = os.path.join(data_output_path, "Union", stat_type)
    if not os.path.exists(union_data_output_dir):
        os.makedirs(union_data_output_dir)
    union_obj = Union_Data(instance_Merged_Data, instance_Exp_Data)

    # Save results to disk
    union_obj.save_union(union_data_output_dir)
    union_obj.save_summary_union(union_data_output_dir)


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
        "diff_ex_dir", type=str, help="parent path of input directory",
    )

    parser.add_argument(
        "stat_type", type=str, help="what type of error correction did you use",
    )

    parser.add_argument(
        "--output_directory",
        type=str,
        help="parent path of output directory",
        default=os.path.join(path_main, "../../../Blueberry_Data/AtBB/data_output"),
    )

    parser.add_argument(
        "--input_directory",
        type=str,
        help="parent path of input directory",
        default=os.path.join(path_main, "../../../Blueberry_Data/AtBB/data_input"),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.homolog_input_file = os.path.abspath(args.homolog_input_file)
    args.diff_ex_dir = os.path.abspath(args.diff_ex_dir)
    args.output_directory = os.path.abspath(args.output_directory)
    args.input_directory = os.path.abspath(args.input_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    logger.info("Starting filtration...")
    process(
        args.syntelog_input_file,
        args.homolog_input_file,
        args.output_directory,
        args.diff_ex_dir,
        args.stat_type,
    )

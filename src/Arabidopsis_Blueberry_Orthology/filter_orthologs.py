#!/usr/bin/env python3

"""
Master code file.
Filter homolog (BLAST) data
Filter syntelog (SynMap) data
Merge both files together
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs

from import_syntelogs import import_syntelogs
from import_homologs import import_homologs
from merge_homo_synt import merge_homo_synt


def process(
    syntelog_input_file,
    homolog_input_file,
    data_output_path,
):
    # Import the synteny data from raw file
    logger.info("Importing syntelogs: %s" % syntelog_input_file)
    syntelogs = import_syntelogs(syntelog_input_file)
    # Save syntelog data to disk
    # MAGIC filename
    file_to_save = os.path.join(data_output_path, "Filtered_Syntelogs.tsv")
    logger.info("Writing syntelog data to disk: %s" % file_to_save)
    syntelogs.to_csv(file_to_save, sep="\t", header=True, index=False)

    # Import the homology (BLAST) data from raw file
    logger.info("Importing BLAST results: %s" % homolog_input_file)
    homologs = import_homologs(homolog_input_file)
    # Save homolog data to disk
    # MAGIC filename
    file_to_save = os.path.join(data_output_path, "Filtered_Homologs.tsv")
    logger.info("Writing homolog data to disk: %s" % file_to_save)
    homologs.to_csv(file_to_save, sep="\t", header=True, index=False)

    # NOTE
    # For both the synteny set and the homology set, if there was a situation
    # where a blueberry gene pointed to multiple non-unique Arabidopsis genes I
    # only took the Blueberry-Arabidopsis pair in which the E-value was
    # greatest. This makes it so that both dataframes are composed of
    # non-unique blueberry genes, Arabidopsis genes may repeat however.

    # Merge the synteny and homology data
    logger.info("Merging the data...")
    merged_all = merge_homo_synt(syntelogs, homologs)
    # Save synteny/homology data to disk
    # MAGIC filename
    file_to_save = os.path.join(data_output_path, "Synteny_Homology_Table.tsv")
    logger.info("Writing merged data to disk: %s" % file_to_save)
    merged_all.to_csv(file_to_save, sep="\t", header=True, index=False)


if __name__ == "__main__":
    """Command line interface to link orthologs together."""

    parser = argparse.ArgumentParser(description="Filter orthologs")
    path_main = os.path.abspath(__file__)
    parser.add_argument(
        "syntelog_input_file", type=str, help="parent path of syntelog file"
    )
    parser.add_argument(
        "homolog_input_file", type=str, help="parent path of homolog file"
    )

    parser.add_argument(
        "output_directory",
        type=str,
        help="parent path of output directory",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.homolog_input_file = os.path.abspath(args.homolog_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    logger.info("Starting filtration...")
    process(
        args.syntelog_input_file,
        args.homolog_input_file,
        args.output_directory,
    )

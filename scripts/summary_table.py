#!/usr/bin/env python3

"""
Create the union of dataframes of differentially expressed genes, Arabidopsis
ortholog, and gene network module identity.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import pandas as pd
import re
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    # TODO get argument for output module information file (blueberry genes)
    parser.add_argument("", type=str, help="")

    # TODO get argument for output Synteny-Homology file/table containing
    # Blueberry:Arabidopsis gene pairs
    parser.add_argument("", type=str, help="")

    # TODO get argument for output differential expression folder
    # NOTE all haplotypes, FDR corrected
    # NOTE 'direction' files will have to be iterated over
    parser.add_argument("", type=str, help="")

    # TODO get argument for output file containing Arabidopsis gene and GO info
    parser.add_argument("", type=str, help="")

    parser.add_argument("output_dir", type=str, help="parent path to output results")
    args = parser.parse_args()

    # TODO redefine all input args as abspath
    # args.go_master_file = os.path.abspath(args.go_master_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE code from here

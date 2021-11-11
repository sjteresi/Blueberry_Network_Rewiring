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

# TODO create function to read genes_and_module_colors file

# TODO create function to read synteny_homology_table file

# TODO create function to iterate through diff_ex_dir and read direction files

# TODO create function to read arabidopsis_go_terms file

# TODO create function to merge all the datasets

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

    parser.add_argument("output_dir", type=str, help="parent path to output results")

    args = parser.parse_args()
    args.genes_and_module_colors = os.path.abspath(args.genes_and_module_colors)
    args.synteny_homology_table = os.path.abspath(args.synteny_homology_table)
    args.diff_ex_dir = os.path.abspath(args.diff_ex_dir)
    args.arabidopsis_go_terms = os.path.abspath(args.arabidopsis_go_terms)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # NOTE code from here

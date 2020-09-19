#!/usr/bin/env python3

"""
Calculates the proportions of each gene belonging to each identification type.
"""

__author__ = "Alder Fulton"

import argparse
import os
import pandas as pd





def process(selection
    """"
    Args:
        selection (str): string of the path to gene annotation file. GFF format.

    Returns:
        None. Prints the proportions of each gene.

    """









if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate gene proportions")

    parser.add_argument("input_file", type=str, help='parent path to the gene file')
    parser.add_argument("selection", type=str, help='all or single haplotype')

    args = parser.parse_args()
    args.input_file = os.path.abspath(args.count_matrix) #Figure this part out.



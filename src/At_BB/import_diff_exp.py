#!/usr/bin/env python3

"""
Import the differential expression data and provide a class for its access
"""

__author__ = "Scott Teresi"

import os
import pandas as pd


class ExpData(object):
    def __init__(self, differential_exp_dir):
        """Initialize
        Args:
            differential_exp_dir (str): String representing the path of the
            directory containing the differentially expressed genes.
            They are outputted from EdgeR.
        """
        self.exp_dict = self.load_diff_exp_sets(differential_exp_dir)

    @staticmethod
    def load_diff_exp_sets(differential_exp_dir):
        """
        Args:
            differential_exp_dir (str): String representing the path of the
            directory containing the differentially expressed genes.
            They are outputted from EdgeR.

        Returns:
            my_file_dict (dictionary of pandaframes):
        """
        my_file_dict = {}
        for entry in os.scandir(differential_exp_dir):
            if entry.path.endswith("Direction.tsv") and entry.is_file():
                my_file = entry.path
                pandas_dataframe = pd.read_csv(my_file, sep="\t", header="infer")

                pandas_dataframe = ExpData.filter_diff_exp(pandas_dataframe)

                my_file_dict[
                    os.path.splitext(os.path.basename(my_file))[0]
                ] = pandas_dataframe
        return my_file_dict

    @staticmethod
    def filter_diff_exp(diff_exp_dataframe):
        """
        Filter the differentially expressed genes, values can be either -1, 0,
        or 1. Non-zero tagged genes were expressed in a positive or negative
        direction and are significant

        Args:
            diff_exp_dataframe (pandas DataFrame): See above
        Returns:
            diff_exp_dataframe (pandas DataFrame): Now with only rows that had
            a -1 or 1"""
        return diff_exp_dataframe.loc[
            diff_exp_dataframe["Direction_Differentially_Regulated"] != 0
        ]

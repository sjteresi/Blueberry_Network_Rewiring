#!/usr/bin/env python3

"""
Import the differential expression data and provide a class for its access
"""

__author__ = "Scott Teresi"

import os
import pandas as pd


class ExpData(object):
    def __init__(self, differential_exp_dir):
        """Initialize"""
        self.exp_dict = self.load_diff_exp_sets(differential_exp_dir)

    @classmethod
    def load_diff_exp_sets(cls, differential_exp_dir):
        """
        Args:
            differential_exp_dir (str): Directory of datafiles

        Returns:
            my_file_dict (dictionary of pandaframes):
        """
        my_file_dict = {}
        for entry in os.scandir(differential_exp_dir):
            if entry.path.endswith("Direction.tsv") and entry.is_file():
                my_file = entry.path
                pandas_dataframe = pd.read_csv(my_file, sep="\t", header="infer")

                pandas_dataframe = cls.filter_diff_exp_p_val(pandas_dataframe)

                my_file_dict[
                    os.path.splitext(os.path.basename(my_file))[0]
                ] = pandas_dataframe
        return my_file_dict

    @staticmethod
    def filter_diff_exp_p_val(diff_exp_dataframe):
        """
        Filter the differential expression data by a P-Value cutoff

        Args:
            differential_exp_dir (str): Directory of datafiles

        Returns:
            diff_exp_dataframe (pandas DataFrame): Now with only rows whose P-Value
            is less than the cutoff
        """
        return diff_exp_dataframe.loc[
            diff_exp_dataframe["Direction_Differentially_Regulated"] != 0
        ]

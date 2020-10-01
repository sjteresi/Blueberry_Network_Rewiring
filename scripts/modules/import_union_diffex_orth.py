#!/usr/bin/env python3

"""
Import differentially expressed genes which were already filtered based on if
they had a matching Arabidopsis gene or not.
"""

__author__ = "Scott Teresi"

import logging
import pandas as pd
import os
import glob


class DiffEx_Orth(object):
    """

    """

    def __init__(self, diffex_orth_dir, logger=None):
        """

        """
        self._logger = logger or logging.getLogger(__name__)
        self.diffex_orth_dir = diffex_orth_dir

        # dirs = os.listdir(diffex_orth_dir)
        # for dirs, root, files in os.walk(diffex_orth_dir):
        # for a_file in files:
        # if a_file.startswith("Union_no_dups"):
        # print(a_file)

    @staticmethod
    def read_diffex_orth(filepath):
        dataframe = pd.read_csv(filepath, header="infer", sep="\t")
        return dataframe

    @property
    def no_duplicates_list(self):
        paths = []
        for dirs, root, files in os.walk(self.diffex_orth_dir):
            for a_file in files:
                if a_file.startswith("Union_no_dups"):
                    paths.append(os.path.join(dirs, a_file))

        return paths

    @property
    def yes_duplicates(self):
        paths = []
        for dirs, root, files in os.walk(self.diffex_orth_dir):
            for a_file in files:
                if a_file.startswith("Union_w_dups"):
                    paths.append(os.path.join(dirs, a_file))

        return paths

    @staticmethod
    def import_data(genes_w_module_groups):
        """

        """
        dataframe = pd.read_csv(
            genes_w_module_groups, sep="\t", header="infer", engine="python",
        )

        return dataframe

    def write_dataframe(self, output_dir, filename):
        """

        """
        self.dataframe.to_csv(
            os.path.join(output_dir, filename), sep="\t", header=True, index=True
        )

    # def __repr__(self):
    # """String representation for developer."""

    # return "Gene_Modules{}".format(self.dataframe.head())

#!/usr/bin/env python3

"""
Wrappers for input data, multiple syntelog pairs.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Scott Teresi"

import logging
import numpy as np
import pandas as pd


class Syntelog_Data(object):
    """

    """

    def __init__(self, syntelog_dataframe, logger=None):
        """Initialize.

        Args:
            syntelog_dataframe (DataFrame): syntelog data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = syntelog_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, blueberry in right-hand column.

        Args:
            filename (str):
        """
        self.to_write = self.dataframe.copy(deep=True)
        # self.to_write.drop(
        # columns=["E_Value", "OrgA_Chromosome", "OrgB_Chromosome", "Diagonal_Score"],
        # inplace=True,
        # )
        self.to_write.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Syntelog_Data{}".format(self.dataframe)

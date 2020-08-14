#!/usr/bin/env python3

"""
Wrappers for ouput data, multiple syntelog and ortholog pairs.
"""

__author__ = "Scott Teresi"

import logging
import numpy as np
import pandas as pd


class Merged_Data(object):
    """

    """

    def __init__(self, merged_dataframe, logger=None):
        """Initialize.

        Args:
            merged_dataframe (DataFrame): data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = merged_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, blueberry in right-hand column.

        Args:
            filename (str):
        """
        self.to_write = self.dataframe.copy(deep=True)
        self.to_write.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Merged_Data{}".format(self.dataframe)

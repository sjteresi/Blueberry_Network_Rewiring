#!/usr/bin/env python3

"""
Wrappers for input data, multiple ortholog pairs.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Scott Teresi"

import logging
import numpy as np
import pandas as pd


class Ortholog_Data(object):
    """

    """

    def __init__(self, ortholog_dataframe, logger=None):
        """Initialize.

        Args:
            ortholog_dataframe (DataFrame): ortholog data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = ortholog_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, blueberry in right-hand column.

        Args:
            filename (str):
        """
        self.dataframe.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Ortholog_Data{}".format(self.dataframe)

#!/usr/bin/env python3

"""
Import gene pairs output, which is dataframe of blueberry and Arabidopsis genes
"""

__author__ = "Scott Teresi"

import logging
import pandas as pd
import os


class Gene_Pairs(object):
    """

    """

    def __init__(self, gene_pairs, logger=None):
        """

        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = self.import_data(gene_pairs)

    @staticmethod
    def import_data(gene_pairs):
        """

        """
        dataframe = pd.read_csv(gene_pairs, sep="\t", header="infer", engine="python",)
        return dataframe

    @property
    def blueberry_gene_list(self):
        """

        """
        return self.dataframe.Blueberry_Gene.tolist()

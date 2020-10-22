#!/usr/bin/env python3

"""
Import WGCNA output, which is dataframe of genes with their module grouping.
"""

__author__ = "Scott Teresi"

import logging
import pandas as pd
import os


class Gene_Modules(object):
    """

    """

    def __init__(self, genes_w_module_groups, logger=None):
        """

        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = self.import_data(genes_w_module_groups)

    @staticmethod
    def import_data(genes_w_module_groups):
        """

        """
        dataframe = pd.read_csv(
            genes_w_module_groups, sep="\t", header="infer", engine="python",
        )

        return dataframe

    @property
    def blueberry_gene_list(self):
        """

        """
        return self.dataframe.Gene_Names.tolist()

    def write_dataframe(self, output_dir, filename):
        """

        """
        self.dataframe.to_csv(
            os.path.join(output_dir, filename), sep="\t", header=True, index=True
        )

    def __repr__(self):
        """String representation for developer."""

        return "Gene_Modules{}".format(self.dataframe.head())

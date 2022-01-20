#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import logging
import os


class Union_Data(object):
    """ """

    def __init__(self, merged_data, exp_data, logger=None):
        """Initialize.

        Args:
            merged_data (object): data frame.
            exp_data (object): data frame.
            output_dir (str):
        """
        self._logger = logger or logging.getLogger(__name__)
        self.merged = merged_data
        self.exp_dict = exp_data.exp_dict
        self.union_w_dups_dict = self.find_union(True)
        self.union_no_dups_dict = self.find_union(False)

    @staticmethod
    def verify_dir(input_filename, output_dir):
        """ """
        the_path = input_filename.replace("_Direction", "")
        if not os.path.exists(os.path.join(output_dir, the_path)):
            os.makedirs(os.path.join(output_dir, the_path))

    def save_union(self, output_dir):
        """

        Args:
            output_dir (str):
        """
        for the_filename, sig_genes in self.union_no_dups_dict.items():
            Union_Data.verify_dir(the_filename, output_dir)
            loc = os.path.join(
                output_dir, the_filename, str("Union_no_dups_" + the_filename + ".tsv")
            )
            self._logger.info("Saving the union data at: %s" % loc)
            sig_genes["Arabidopsis_Gene"].to_csv(
                loc,
                sep="\t",
                header=True,
                index=False,
            )

        for the_filename, sig_genes in self.union_w_dups_dict.items():
            Union_Data.verify_dir(the_filename, output_dir)
            loc = os.path.join(
                output_dir, the_filename, str("Union_w_dups_" + the_filename + ".tsv")
            )
            self._logger.info("Saving the union data at: %s" % loc)
            sig_genes["Arabidopsis_Gene"].to_csv(
                loc,
                sep="\t",
                header=True,
                index=False,
            )

    def save_summary_union(self, output_dir):
        """
        Args:
            output_dir (str):
        """
        for the_filename, sig_genes in self.union_w_dups_dict.items():
            the_filename = the_filename.replace("_Direction", "")
            total_arabidopsis_genes = sig_genes.shape[0]
            unique_arabidopsis_genes = len(sig_genes.Arabidopsis_Gene.unique())
            try:
                ratio = unique_arabidopsis_genes / total_arabidopsis_genes
            except ZeroDivisionError:
                ratio = 0
            summary_info = {
                "Unique_Arabidopsis": unique_arabidopsis_genes,
                "Total_Arabidopsis": total_arabidopsis_genes,
                "Ratio_U/T": ratio,
            }
            summary_frame = pd.DataFrame.from_dict(
                summary_info, columns=["Value"], orient="index"
            )

            the_filename = the_filename.replace("_Direction", "")
            loc = os.path.join(
                output_dir, the_filename, str("Union_Summary_" + the_filename + ".tsv")
            )
            self._logger.info("Saving the summary data at: %s" % loc)
            summary_frame.to_csv(
                os.path.join(
                    output_dir,
                    the_filename,
                    str("Union_Summary_" + the_filename + ".tsv"),
                ),
                sep="\t",
                header=True,
            )

    def find_union(self, dups):
        """
        Args:
            dups (bool): Keep duplicates if true
        """
        the_union_dict = {}
        for the_filename, sig_genes in self.exp_dict.items():
            the_filename = the_filename.replace("_Direction", "")
            diffexp_genes_list = sig_genes.Gene_Name.tolist()

            AT_union = self.merged[
                self.merged["Blueberry_Gene"].isin(diffexp_genes_list)
            ]
            # NOTE sort by E_Value so that when we do dorp duplicates, we can
            # keep the first, which is the best score (most smallest).
            AT_union = AT_union.sort_values(by=["E_Value"])
            if dups:
                the_union_dict[the_filename] = AT_union
            else:
                the_union_dict[the_filename] = AT_union.drop_duplicates(
                    subset=["Arabidopsis_Gene"], keep="first"
                )
        return the_union_dict

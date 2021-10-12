#!/usr/bin/env python3

"""
Functions to merge the homology (BLAST) and synteny (SynMap) output that has
been pre-filtered. Helper functions identify genes that are in the homolog data
(more complete)
"""

__author__ = "Scott Teresi"

import pandas as pd


def merge_homo_synt(syntelogs, homologs):
    """
    Only add a blueberry gene from the BLAST set (wrapped homologs) if that
    gene is not currently in the synteny set. This is because we give priorty
    to results determined through synteny.

    Args:
        syntelogs (Syntelog_Data object): wrapped data

        homologs (Homolog_Data object ): wrapped data

    Returns:
        merged_all (pandas.core.Data.Frame): Pandas dataframe that is the
        amalgamation of the syntelog and homolog data.
    """
    syntelog_blueberry_list = syntelogs.Blueberry.to_list()
    boolean_mask = homologs.Blueberry.isin(syntelog_blueberry_list)
    missing_links = homologs[~boolean_mask]

    merged_all = pd.concat([syntelogs, missing_links], axis=0, join="outer")
    merged_all.rename(
        columns={"Arabidopsis": "Arabidopsis_Gene", "Blueberry": "Blueberry_Gene"},
        inplace=True,
    )
    return merged_all

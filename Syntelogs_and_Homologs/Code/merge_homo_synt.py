#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd


def merge_homo_synt(wrapped_syntelogs, wrapped_homologs):
    """
    Only add a blueberry gene from the BLAST set (wrapped homologs) if that
    gene is not currently in the synteny set.
    """
    missing_links = _identify_missing_links(wrapped_syntelogs, wrapped_homologs)
    merged_all = pd.concat(
        [wrapped_syntelogs.dataframe, missing_links], axis=0, join="outer"
    )
    merged_all.rename(
        columns={"Arabidopsis": "Arabidopsis_Gene", "Blueberry": "Blueberry_Gene"},
        inplace=True,
    )
    return merged_all


def _identify_missing_links(wrapped_syntelogs, wrapped_homologs):
    """
    Create and use a boolean mask to identify the blueberry genes in the
    homolog data that are missing from the syntelog data.

    """
    syntelog_blueberry_list = wrapped_syntelogs.dataframe.Blueberry.to_list()
    boolean_mask = wrapped_homologs.dataframe.Blueberry.isin(syntelog_blueberry_list)
    missing_links = wrapped_homologs.dataframe[~boolean_mask]
    return missing_links

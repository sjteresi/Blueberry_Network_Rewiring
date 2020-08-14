import pandas as pd


def merge_ortho_synt(wrapped_syntelogs, wrapped_orthologs):
    """
    Only add a blueberry gene from the BLAST set (wrapped orthologs) if that
    gene is not currently in the synteny set.
    """
    missing_links = _identify_missing_links(wrapped_syntelogs, wrapped_orthologs)
    merged_all = pd.concat(
        [wrapped_syntelogs.dataframe, missing_links], axis=0, join="outer"
    )
    merged_all.rename(
        columns={"Arabidopsis": "Arabidopsis_Gene", "Blueberry": "Blueberry_Gene"},
        inplace=True,
    )
    return merged_all


def _identify_missing_links(wrapped_syntelogs, wrapped_orthologs):
    """

    """
    syntelog_blueberry_list = wrapped_syntelogs.dataframe.Blueberry.to_list()
    boolean_mask = wrapped_orthologs.dataframe.Blueberry.isin(syntelog_blueberry_list)
    missing_links = wrapped_orthologs.dataframe[~boolean_mask]
    return missing_links

"""
Read files and directories. Perform minimal re-formatting on the files.
If this should be done, it will be done in the main analysis scripts.
"""

__author__ = "Scott Teresi"

import argparse
import coloredlogs
import logging
import numpy as np
import pandas as pd
import os


def read_GO_dir(input_directory):
    """
    Reads an input directory and returns a list of file paths that contain
    the magic substring 'Summary'

    Args:
        input_directory (str): Path to directory that contains files with our
        magic substring. This directory is parsed to find filenames that match
        the substring.

    Returns:
        go_files (list): List of absolute file paths
    """
    go_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("Summary" in f)
    ]  # MAGIC substring for filename recognition
    return go_files

def read_protein_table(filepath)
    """
    Read the table of Arabidopsis genes and protein ID/names that was
    previously created.

    Args:
        filepath (str): Path to the table

    Returns:
        protein_table (pandas dataframe): columns=[Arabidopsis_Gene,
        Protein_ID, Protein_Name]
    """
    protein_table = pd.read_csv(filepath, sep='\t', header='infer')
    return protein_table


def read_gene_modules_table(filepath):
    """
    Read the table of blueberry genes and module identities that was
    previously created from WGCNA

    Args:
        filepath (str): Path to the table

    Returns:
        gene_modules_table (pandas dataframe): columns=[Blueberry_Gene,
        Module_Color]
    """
    gene_modules_table = pd.read_csv(
        filepath,
        sep="\t",
        header="infer",
    )
    gene_modules_table.rename(
        columns={"Gene_Names": "Blueberry_Gene", "moduleColor": "Module_Color"},
        inplace=True,
    )

    return gene_modules_table


def read_DEG_dir(input_directory):
    """
    Reads an input directory and returns a list of file paths that contain
    the magic substring 'Direction.tsv'.

    Args:
        input_directory (str): Path to directory that contains files with our
        magic substring. This directory is parsed to find filenames that match
        the substring.

    Returns:
        deg_files (list): List of absolute file paths
    """
    deg_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if (os.path.isfile(os.path.join(input_directory, f))) and ("Direction" in f)
    ]  # MAGIC substring for filename recognition

    if len(deg_files) != 28:  # MAGIC number of files, anymore or less will
        # result in error
        raise ValueError(
            """More files were tagged than expected, check to see
                         if files are being tagged appropriately."""
        )

    return deg_files


def read_all_DEG_direction_files(deg_files):
    """
    Reads a list of file paths and reads them as pandas files. Creates a list
    of pandas dataframes. Magic column names are read from the file

    Args:
        deg_files (str): Output from read_log_2fc_dir function. List
        of file paths.

    Returns:
        list_of_read_deg_files (list): List of pandas dataframes
         TODO fill in to be more descriptive for the structure
    """
    list_of_read_deg_files = []
    for deg_path in deg_files:
        individual_deg_file = pd.read_csv(deg_path, sep="\t", header="infer")

        if len(os.path.splitext(os.path.basename(deg_path))) > 2:
            # MAGIC remove filename extension, does not work if multiple periods
            # are in filename.
            raise ValueError(
                """There were multiple periods in your file path,
                             this will cause an error."""
            )
        analysis_context = os.path.splitext(os.path.basename(deg_path))[0]
        analysis_context = analysis_context.rstrip("_Direction")  # NB strip
        # the string from the end so the groupby unstack command in
        # convert_direction_integer_to_string gives us a more appropriate
        # column name

        individual_deg_file.rename(
            columns={
                "Direction_Differentially_Regulated": analysis_context,
                "Gene_Name": "Blueberry_Gene",
            },
            inplace=True,
        )  # MAGIC column name
        individual_deg_file = convert_direction_integer_to_string(
            individual_deg_file, analysis_context
        )
        list_of_read_deg_files.append(individual_deg_file)

    return list_of_read_deg_files


def convert_direction_integer_to_string(deg_panda_file, analysis_context):
    """
    Take the output from EdgeR (the direction files), where a gene can either
    be 0, 1, or -1 and convert those integer files to strings signifying the
    direction of regulation (i.e up, down, or no change)

    Args:
        deg_panda_file (pandas.DataFrame):
            Index:
                RangeIndex
            Columns:
                Blueberry_Gene: Object, dtype: object
                analysis_context: Integer, dtype: int64
            Shape:
                (93583, 2)

        analysis_context (string):
            String representing the column name and comparison context
            of the EdgeR data.

    Returns: deg_panda_file (pandas.DataFrame):
        Index:
            RangeIndex
        Columns:
            Blueberry_Gene: Object, dtype: object
            analysis_context + _'Direction': Object, dtype: object
        Shape:
            (93583, 2)
    """
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == -1, analysis_context + "_Direction"
    ] = "Down"
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == 1, analysis_context + "_Direction"
    ] = "Up"
    deg_panda_file.loc[
        deg_panda_file[analysis_context] == 0, analysis_context + "_Direction"
    ] = "No_Change"
    deg_panda_file.drop(columns=analysis_context, inplace=True)

    return deg_panda_file


def read_GO_ID_w_term(input_file):
    """
    Read a GO table from file
    """
    data = pd.read_csv(
        input_file,
        header="infer",
        sep="\t",
        dtype={"GO_Term": "object", "GO_ID": "object"},
    )
    return data


def read_synteny_homology_table(filepath):
    """
    Read the synteny/homology table of blueberry and arabidopsis genes that was
    previously created
    Args:
        filepath (str): Path to the table

    Returns:
        synteny_homology_table (pandas dataframe): columns = [
        Arabidopsis_Gene, Blueberry_Gene, E_Value, Point_of_Origin]
    """
    synteny_homology_table = pd.read_csv(
        filepath,
        sep="\t",
        header="infer",
    )
    return synteny_homology_table


def read_interesting_modules(filepath):
    """
    TODO
    """
    colnames = ["Modules"]
    interesting_modules = pd.read_csv(filepath, names=colnames, header=None, sep="\t")
    return interesting_modules


def read_interesting_GO_terms(filepath):
    """
    TODO
    """
    colnames = ["GO_ID"]
    interesting_go_terms = pd.read_csv(filepath, names=colnames, header=None, sep="\t")
    return interesting_go_terms


def read_FPKM_or_TPM(filepath):
    data = pd.read_csv(filepath, header="infer", sep="\t")
    data.rename(columns={"Gene_Name": "Blueberry_Gene"}, inplace=True)
    return data

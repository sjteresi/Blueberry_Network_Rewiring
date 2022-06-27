#!/usr/bin/env python3

"""
Creates a histogram of the percent expression of modules in each comparison.
"""

__author__ = "Alder Fulton"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "file_path",
        type=str,
        help="""
                        The path to filtered log 2fc data file""",
    )
    parser.add_argument("output_path", type=str, help="The output directory")
    args = parser.parse_args()
    file_path = args.file_path
    output = args.output_path

    data = pd.read_csv(file_path, delimiter="\t")

    # NB create multiplot obj
    fig, axes = plt.subplots(5, 6, figsize=(30, 30 * 5 / 6))
    axes = np.array(axes).ravel()
    fig.patch.set_facecolor("white")

    # NOTE or TODO this creates a square figure of 30 subplots but 2 remain
    # unfilled. Should this be edited? Or because it will be supplemental does
    # it matter?
    i = 0
    for col in data.columns[1:29]:  # Magic, use columns of filtered data. Make
        # subplot for each column
        axes[i].title.set_text(col)
        axes[i].hist(data[col])
        i += 1

    # TODO this could probably be explained better. Reword?
    # NB likely a supplemental figure, maybe get Pat's opinion.
    fig.suptitle(
        """Histogram of the number of modules and the percentage of genes in
        each log2FC expression context that overlap with each module""",
        fontsize=60,
    )
    # TODO double-check?
    fig.supxlabel("Percent of genes recovered by context per module", fontsize=50)
    fig.supylabel("Number of Modules", fontsize=50)

    plt.savefig(output)

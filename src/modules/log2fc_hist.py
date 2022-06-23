#!/usr/bin/env python3

"""
Creates a histogram of the percent exrpression of modules in each comparison.
"""

__author__ = "Alder Fulton"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("file_path", type=str, help="The path to the file")
    parser.add_argument("output_path", type=str, help="The path to save the file")
    args = parser.parse_args()
    file_path = args.file_path
    output = args.output_path

    data = pd.read_csv(file_path, delimiter="\t")

    fig, axes = plt.subplots(5, 6, figsize=(30, 30 * 5 / 6))
    axes = np.array(axes).ravel()
    fig.patch.set_facecolor("white")

    i = 0
    for col in data.columns[1:29]:  # Magic
        axes[i].title.set_text(col)
        axes[i].hist(data[col])
        i += 1
    fig.suptitle("Distribution of connection strengths\n between modules and comparisons", fontsize=80)
    fig.supxlabel("Percent Expression", fontsize=50)
    fig.supylabel("Count", fontsize=50)

    plt.savefig(output)

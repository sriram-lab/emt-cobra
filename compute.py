"""
"""

import pandas as pd
import numpy as np


def log2FoldChange(df):
    """
    """
    control = df[["Intensity CTRL1", "Intensity CTRL2"]].mean(axis=1)
    treatment = df[["Intensity TREATED1", "Intensity TREATED2"]].mean(axis=1)

    foldChange = treatment.div(control)
    maxFC = foldChange.loc[[foldChange != np.inf]].max()
    foldChange = foldChange.replace(np.inf, maxFC, inplace=True)

    log2FoldChange = foldChange.apply(np.log2)

    df["log2 FC"] = log2FoldChange
    return df


def computeZScore(df):
    """
    """
    mean = df["log2 FC"].mean()
    stddev = df["log2 FC"].std()
    diff = df["log2 FC"] - mean
    df["Zscore"] = diff.div(stddev)
    return df

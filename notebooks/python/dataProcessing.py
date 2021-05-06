"""
@author: Scott Campit
"""

import pandas as pd
import numpy as np
import scipy.stats

def log2FC(data1, data2):
    """
    """

    foldChange = data1.div(data2)
    log2FoldChange = foldChange.apply(np.log2)

    return log2FoldChange

def oneSampleZTest(sample):
    """
    """

    mu = sample.mean()
    var = sample.std()
    diff = (sample - mu)
    zscore = diff.div(var)
    pvalue = scipy.stats.norm.sf(abs(zscore))*2

    return zscore, pvalue

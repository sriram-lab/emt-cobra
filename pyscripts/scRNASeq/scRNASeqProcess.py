"""
splitting scRNA data into 4 separate files

"""

import sys
import pandas as pd
import numpy as np


def split_dataset(file, delimit='\t'):
    df = pd.read_csv(file, delimiter=delimit)
    conditions = ["day8_rep1", "day8_rep2", "day10_rep1", "day10_rep2"]

    for cond in conditions:
        tmp = df[df["day_replicate"] == cond]
        tmp.to_csv(cond + ".csv", index=False)


file = sys.argv[1]
split_dataset(file)

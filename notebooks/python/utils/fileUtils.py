"""
"""

import sys
import pandas as pd

def tsv_to_csv(file):
    """
    """
    name = file.split(".")[0]
    df = pd.read_csv(file, delimiter='\t')
    df.to_csv(name+'.csv')

if __name__ == "__main__":
    tsv_to_csv(sys.argv[1])

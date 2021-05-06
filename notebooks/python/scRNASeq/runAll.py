"""
"""

import pandas as pd
import mapping
import scRNASeqProcess

file = "./../../datasets/scRNASeq/GSE114397_HMLE_TGFb.tsv"
df = pd.read_csv(file, delimiter='\t')
df = df.transpose()
print(df)
# Get mapped gene dataset

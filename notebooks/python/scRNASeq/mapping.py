"""
mapping data
"""

import sys
import pandas as pd
import mygene


def getGeneSymbols(OriginalDataSet, OriginalGeneIDs):
    mg = mygene.MyGeneInfo()
    geneName = list(OriginalGeneIDs.values.flatten())

    geneDB = mg.querymany(geneName,
                          scopes=['ensemble', 'ensemble.transcript'],
                          fields='symbol', species='Human', returnall=True, as_dataframe=False)
    print(geneDB.values())
    #OriginalDataSet = OriginalDataSet.rename(list(geneDB.index), axis=1)
    #print(OriginalDataSet)
    #return OriginalDataSet


tgfdata = pd.read_csv(
          './../../datasets/scRNASeq/GSE114397_HMLE_TGFb.tsv', delimiter='\t')
originalIDs = tgfdata.columns

MappedSet = getGeneSymbols(tgfdata, originalIDs)

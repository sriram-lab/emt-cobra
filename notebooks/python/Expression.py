"""
@author: Scott Campit
"""

import mygene
import pandas as pd
import numpy
import GEOparse
import compute


def mapGenesToData(file):
    """
    """
    gse = GEOparse.get_GEO(file, destdir='./')
    for _, gsm in gse.gsms.items():
        accession = gsm.table.set_index('ID_REF')
    data = gse.pivot_samples('VALUE')
    data = pd.merge(data, accession,
                    how='inner',
                    left_index=True, right_index=True)

    data = data.filter(regex='(GSM)')
    for _, gpl in gse.gpls.items():
        GPLTable = gpl.table.set_index("ID")

    MappedSet = pd.merge(data, GPLTable,
                         how='inner',
                         left_index=True, right_index=True)

    if "GB_ACC" in MappedSet.columns:
        geneName = MappedSet.GB_ACC
        MappedSet = MappedSet.set_index('GB_ACC')
    else:
        pass

    return MappedSet, geneName


def getGeneSymbols(OriginalDataSet, OriginalGeneIDs, Query):
    mg = mygene.MyGeneInfo()
    geneName = list(OriginalGeneIDs.values.flatten())

    if Query == 'DNA_Array':
        geneDB = mg.querymany(geneName,
                              scopes=['accession.rna', 'refseq', 'reporter',
                                      'entrezgene', 'hprd'],
                              fields='symbol', species='Human', as_dataframe=True)

    elif Query == 'RNA_Seq':
        geneDB = mg.querymany(geneName,
                              scopes=['ensemble', 'ensemble.transcript'],
                              fields='symbol', species='Human', as_dataframe=True)

    MappedSet = pd.merge(OriginalDataSet, geneDB,
                         how='inner',
                         left_index=True, right_index=True)
    MappedSet = MappedSet.filter(
            regex=r'(GSM|symbol)')
    MappedSet = MappedSet.set_index('symbol')

    return MappedSet


def mapGSMtoLabel(ExpData):
    map = pd.read_csv('./../datasets/GSM_map.csv')
    newNames = list(map['Type'].values)
    ExpData.columns = newNames

    return ExpData


#file = r'GSE17708'
#MappedSet, geneName = mapGenesToData(file)
#finalExpData = getGeneSymbols(
#        MappedSet, geneName)
#finalExpData.to_csv('./GSE17708mapped.txt', index=True)
#data = pd.read_csv('GSE17708mapped.txt', index_col='symbol')
#mappedData = mapGSMtoLabel(data)
#mappedData.to_excel('./GSE17708mapped.xlsx', index=True)
df = pd.read_excel('./GSE17708mapped.xlsx', index_col='symbol')
geneNames = list(df.index.values)
conditionNames = list(df.columns.values)
data = df.values
RcoefMat, PvalMat = compute.computeColumnZ(data)
print(RcoefMat)
print(PvalMat)

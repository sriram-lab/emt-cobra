"""

"""
import pandas as pd
import numpy as np


def destack(df, column, delimiter=";"):
    vect = df[column].str.split(delimiter).apply(pd.Series, 1).stack()
    vect.index = vect.index.droplevel(-1)
    vect.name = column
    del df[column]
    df = df.join(vect)
    df = df.drop_duplicates(keep='first')

    return df


df = pd.read_excel('Proteomics2010.xlsx')
newDF = destack(df, "Gene Names", delimiter=';')
newDF = newDF[newDF.Contaminant != '+']
newDF = newDF[newDF.Reverse != '+']
newDF = newDF[["Gene Names", "REFSEQ", "Intensity CTRL1", "Intensity CTRL2",
               "Intensity TREATED1", "Intensity TREATED2", "Reverse", "Contaminant"]]
newDF.to_excel('DestackProteomics2010.xlsx')

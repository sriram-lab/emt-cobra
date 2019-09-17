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

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = lst_cols.split(';')
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res


df = pd.read_excel('Proteomics2010.xlsx')
newDF = explode(df, ["Gene Names", "REFSEQ"], fill_value='')
newDF = newDF[newDF.Contaminant != '+']
newDF = newDF[newDF.Reverse != '+']
newDF = newDF[["Gene Names", "REFSEQ", "Intensity CTRL1", "Intensity CTRL2",
               "Intensity TREATED1", "Intensity TREATED2", "Reverse", "Contaminant"]]
newDF.to_excel('DestackProteomics2010.xlsx')

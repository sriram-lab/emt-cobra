"""
make_table formats data into plotly tables for visualization

@author: Scott Campit
"""
import dash_table
import dash_html_components as html
import dash_core_components as dcc
import dash
import chart_studio.plotly as py

import pandas as pd
import numpy as np
import xlrd

import plotly.graph_objs as go
from plotly.subplots import make_subplots
import chart_studio
import plotly.graph_objs as go


def makeTableStructure(data):
    """
    """
    head = []
    df_values = []
    data_length = data.shape[1]
    data = data.sort_values(by='Zscore')
    head = ["Reaction Name", "Metabolic Pathway",
            "A459 Reaction<br>Flux", "P-value", "Z-score"]

    df_values = [data["Reaction"], data["Metabolic_Pathway"], data["Grate"],
                 data["pValue"], data["Zscore"]]

    table = go.Table(
        header=dict(
            values=head,
            font=dict(
                size=10)
            ),
        cells=dict(
            values=df_values
            )
        )
    return table


def make_table(data):
    """
    """
    posTscore = data.loc[(data['pValue'] <= 0.05)]
    negTscore = data.loc[(data['pValue'] <= 0.05)]

    table = pd.concat([posTscore, negTscore], axis=0, sort=['Zscore'])
    tableStructure = makeTableStructure(table)

    return tableStructure

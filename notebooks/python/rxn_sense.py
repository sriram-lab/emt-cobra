"""
"""
import pandas as pd
import numpy as np
import json

from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

import flask


# Intialize the Flask/Dash application
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = 'metabolomics-pipeline'

# Read in data
df = pd.read_json('./../python/me1_filter.json')
df['day_replicate'] = df['Day'].astype(str)+'_'+df['Replicate'].astype(str)

# Create dataframes to parse data into three heat maps
panc1 = df.loc[(df["Cell Line"] == "PANC1")]
panc1 = panc1.sort_values(by='Value', ascending=False)
tu8902 = df.loc[(df["Cell Line"] == "TU8902")]
tu8902 = tu8902.sort_values(by='Value', ascending=False)
tu8988t = df.loc[(df["Cell Line"] == "TU8988T")]
tu8988t = tu8988t.sort_values(by='Value', ascending=False)

# get values that will be used
knockout = df["Knockout"].unique()
mode = df["Mode"].unique()
type = df["Type"].unique()

"""
Functions
@author: Scott Campit
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

df = pd.DataFrame()
colormap = []

def widget(data=df):
    """
    """
    # get values that will be used
    knockout = df["Knockout"].unique()
    mode = df["Mode"].unique()
    type = df["Type"].unique()

    # Create dynamic parts that will allow client to interact with data
    widgets = dbc.Container(
        [
            # Slider
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            dcc.Dropdown(
                                id='knockout-type',
                                options=[{'label': i, 'value': i}
                                         for i in knockout],
                                value='ME1'
                                )
                            ),
                        width={"size": 3.0,
                               "offset": 0.3}
                        ),
                    dbc.Col(
                        html.Div(
                            dcc.Dropdown(
                                id='mode-type',
                                options=[{'label': i, 'value': i}
                                         for i in mode],
                                value='(+)'
                                )
                            ),
                        width={"size": 3.0,
                               "offset": 0.3}
                        ),
                    dbc.Col(
                        html.Div(
                            dcc.Dropdown(
                                id='type',
                                options=[{'label': i, 'value': i}
                                         for i in type],
                                value='Intracellular'
                                )
                            ),
                        width={"size": 3.0,
                               "offset": 0.3}
                        )
                    ]
                )
            ]
        )

    return widgets


def make_struct(hm_id='type-heatmap', data=df, nam='type', cmap=colormap):
    """
    """
    container = html.Div(
        dcc.Graph(
            id=hm_id,
            figure={
                'data': [(
                    go.Heatmap(
                        x=data['Gene'],
                        y=data['Feature'],
                        z=data['Value'],
                        name=nam,
                        colorscale=cmap)
                    )],
                'layout': go.Layout(
                    title=go.layout.Title(
                        text=('<b>' + data['Cell Line'].iloc[0]+'</b>'),
                        xanchor='right',
                        yanchor='bottom',
                        x=0,
                        font=dict(
                            family='Arial',
                            size=16,
                            color='black'
                            )
                        ),
                    autosize=False,
                    yaxis=dict(
                        automargin=True,
                        tickfont=dict(
                            family='Arial, sans-serif',
                            size=14,
                            color='black'
                            )
                        )
                    )
                },
            config={
                'displayModeBar': False
                }
            )
        )
    return container

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

import make_tables
import dataProcessing


def custom_hover(data):
    """
    """
    hover = []
    number_of_objs = data.shape[0]
    dataType = ''

    if number_of_objs < 3000:
        dataType = 'Genes'
        for _, column in data.iterrows():
            hover_data = 'Gene: '+'{}'.format(column['Gene_symbol'])+'<br>' \
                + 'P-value: '+'{:.5f}'.format(column['pValue'])+'<br>' \
                + 'Z-score: '+'{:.2f}'.format(column['Zscore'])
            hover.append(hover_data)

    else:
        dataType = 'Reactions'
        for _, column in data.iterrows():
            hover_data = 'Reaction: '+'{}'.format(column['Reaction'])+'<br>' \
                + 'Metabolic Pathway: '+'{}'.format(column['Metabolic_Pathway'])+'<br>' \
                + 'P-value: '+'{:.5f}'.format(column['pValue'])+'<br>' \
                + 'Z-score: '+'{:.2f}'.format(column['Zscore'])
            hover.append(hover_data)

    return hover, dataType


def makeScatterPlots(data, hover):
    """
    """
    scatter_plot = go.Box(
        x=data['Zscore'].values,
        marker_size=7,
        jitter=0.5,
        hoverinfo='text',
        text=hover,
        marker=dict(
          line=dict(
                      color="#000000",
                      width=0.4
                )
            )
        )
    return scatter_plot


def makeHistogram(data):
    """
    """
    histogram = go.Histogram(
        x=data['Zscore'].values,
        marker_color='#003D73',
        name='',
        marker={"line": {"color": "#000000", "width": 0.4}}
        )

    return histogram


def makeDifferentialSensitivityPlotVars(data):
    """
    """
    hover, dataType = custom_hover(data)
    scatterplot = makeScatterPlots(data, hover)
    histogram = makeHistogram(data)
    tableStructure = make_tables.make_table(data)

    fig = make_subplots(rows=3, cols=1,
                        shared_xaxes=True, shared_yaxes=False,
                        vertical_spacing=0.05,
                        specs=[
                            [{"type": "table"}],
                            [{"type": "scatter"}],
                            [{"type": "histogram"}],
                            ])

    fig.add_trace(tableStructure, row=1, col=1)
    fig.add_trace(histogram, row=2, col=1)
    fig.add_trace(scatterplot, row=3, col=1)
    return fig, dataType


def trace_asethetics(fig):
    """
    """
    fig.update_traces(opacity=0.9, row=2, col=1)
    fig.update_traces(opacity=0.9, row=3, col=1)
    return fig


def figure_labels(fig, dataType):
    """
    """
    if dataType is 'Genes':
        textTitle = "Number of genes: 1488"
    else:
        textTitle = "Number of reactions: 3743"
    fig.update_yaxes(title_text=textTitle, type="log",
                     tickvals=[1, 10, 100, 1000, 10000], row=2, col=1)
    fig.update_yaxes(showticklabels=False, row=3, col=1)
    fig.update_xaxes(title_text="Z-Score", showgrid=False, row=3, col=1)

    return fig


def make_final_figure(data, fileName):
    """
    """
    fig_struct, dataType = makeDifferentialSensitivityPlotVars(data)
    modified_fig = trace_asethetics(fig_struct)
    final_fig = figure_labels(modified_fig, dataType)
    final_fig.update_layout(title=fileName, showlegend=False,
                            height=1000, width=1750)
    return final_fig


def makeHeatmapDF(cellLine, condition):
    control = pd.read_excel(lower(cellLine)+'.xlsx',
                            sheet_name=cellLine+" - dox")
    treatment = pd.read_excel(lower(cellLine)+'.xlsx',
                              sheet_name=cellLine+" "+condition)

    log2FoldChange = dataProcessing.log2FC(treatment, control)

    return log2FoldChange


def plotHeatmap(log2FoldChange):
    None

def makeVolcanoPlot(control, treatment):
    log2FoldChange = dataProcessing.log2FC(treatment, control)
    zscore, pvalue = dataProcessing.oneSampleZTest(log2FoldChange)
    adjPValue = -1*pvalue.apply(np.log10)

    return zscore, adjPValue

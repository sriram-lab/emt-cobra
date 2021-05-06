"""
differential_sensitivity_plots creates plotly figures for each of the treatment conditions.

@author: Scott Campit
"""

import color_wheels
import make_tables
import make_plots

import xlrd
import numpy as np
import pandas as pd

import chart_studio.plotly as py
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import chart_studio
chart_studio.tools.set_credentials_file(
    username='ScottCampit', api_key='ZwIWubUD2bgOSoYzECmG')


data = pd.read_excel('diffrxnsense1.xlsx')

df = data
cancerName = "A459"
fig = make_plots.make_final_figure(df, cancerName)
#fig.show()
py.plot(fig, file_names=cancerName)

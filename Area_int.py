from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px


def plot_areas(N):

    a=1

    ico2 = Polyhedron.icosahedron(fold = 2)
    fold_2 =np.empty(N)

    ico3 = Polyhedron.icosahedron(fold = 3)
    fold_3 =np.empty(N)
    
    ico5 = Polyhedron.icosahedron(fold = 5)
    fold_5 =np.empty(N) 

    interval = np.linspace(0.01,2*a*np.sin(2*np.pi/5),N)

    for i in range(len(interval)):
        fold_2[i] = ico2.area_ico(interval[i],a)
        fold_3[i] = ico3.area_ico(interval[i],a)
        fold_5[i] = ico5.area_ico(interval[i],a)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=interval, y=fold_5, name="5-fold"))
    fig.add_trace(go.Scatter(x=interval, y=fold_3, name="3-fold"))
    fig.add_trace(go.Scatter(x=interval, y=fold_2, name="2-fold"))

    fig.update_layout(
            title='Interacting area of icosahedron vs cutoff height',
            xaxis_title= 'hc (a)',
            yaxis_title="Area (a^2)",
            legend_title="Symmetry:",
            font=dict(family="Courier New, monospace", size=18, color="RebeccaPurple"))
    
    fig.show()

#plot_areas(1000)



























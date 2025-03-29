import plotly.graph_objects as go
import numpy as np
np.random.seed(1)

N = 70
x = np.array([0,1,0,1])
y = np.array([0,0,1,1])
z = np.array([0,0,0,0])
fig = go.Figure(data=[go.Mesh3d(x=x,
                   y=y,
                   z=z,
                   opacity=0.5,
                   color='rgba(244,22,100,0.6)'
                  )])

fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[-100,100],),
                     yaxis = dict(nticks=4, range=[-50,100],),
                     zaxis = dict(nticks=4, range=[-100,100],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))

fig.show()
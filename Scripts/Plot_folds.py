import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import os

if not os.path.exists("images"):
    os.mkdir("images")


def funcM(x,y,hm,k):
    return hm *( np.sin(k*x - 3*np.pi/(2*k) )+np.sin(k*y - 3*np.pi/(2*k)) +2)


# Función que crea el mesh de la molécula
def SupM(funcM,R,M,k,hm):

    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)
    X,Y = np.meshgrid(x,y)
    
    Zm = funcM(X,Y,hm,k)
    return Zm

R = np.sin(2*np.pi/5)

datos = np.loadtxt('icosahedron.txt')
x = datos[:,0]
y = datos[:,1]      
z = datos[:,2] + .1

datos_edges = np.loadtxt('edges.txt')
xe = datos_edges[:,0]
ye = datos_edges[:,1]
ze = datos_edges[:,2] + .1

min_x = np.min(x)
max_x = np.max(x)
min_y = np.min(y)
max_y = np.max(y)
min_z = np.min(z)
max_z = np.max(z)

x_p = np.linspace(2*min_x,2*max_x,100)
y_p = np.linspace(2*min_y,2*max_y,100)
z_p = SupM(funcM,R,100,5,0)


fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(size=1)),go.Surface(x=x_p, y =y_p, z=z_p, showscale=False),
    go.Scatter3d(x=xe, y=ye, z=ze, mode='markers', marker=dict(size=3))])
fig.update_layout(
    title="2-fold symmetry with plane surface",
    font=dict(
        family="Courier New, monospace",
        size=14
    )
)

fig.show()
fig.write_html("/home/carlos/Uni/TFG/Jupyter/2_fold_plane.html")














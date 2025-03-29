from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt
#import plotly.graph_objects as go
#import plotly.express as px


a = 1
N = 100
fold = 5 
ico = Polyhedron.icosahedron(fold = fold)
#ico.plot(show_labels=True)


h_5 = 2*a*np.sin(2*np.pi/5)
h_3 = 2*np.sqrt(3)/12*(3+np.sqrt(5))*a
h_2 = 2*a*np.cos(np.pi/5)


areas = np.empty(N)
hc = a*np.sqrt(1-1/(4*np.sin(np.pi/5)*np.sin(np.pi/5)))
interval = np.linspace(0.1,h_5+0.01,N)
const = np.ones(len(interval))*20*np.sqrt(3)/4
for i in range(len(interval)):
    areas[i] = ico.area_ico(interval[i],a)

plt.plot(interval,areas,'.')
plt.plot(interval,const)
plt.xlabel('hc')
plt.ylabel('area')
plt.show()






































from shapely.geometry import Polygon, Point
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import matplotlib.pyplot as plt

# Defino los parámetros que voy a usar 
a = 1                                                               # Longitud de la arista del icosaedro
F = 2                                                               # Tipo de simetría (X-fold)
R = a*np.sin(2*np.pi/5)                                             # 
N = 50                                                              # Número de iteraciones sobre los parámetros estudiados
M = 100                                                              # Número de puntos creados sobre las caras
hc = a*np.sqrt(1-1/(4*np.sin(np.pi/5)*np.sin(np.pi/5))) - 0.00001*a # Altura de corte para los cálculos (calculada desde 
                                                                    # el punto más bajo del icosahedro)

# Inicializo el volumen (comienza en un valor muy elevado para asegurar que no hay problemas más adelante)
V = 3000*a

# Creo el icosaedro 
p = Polyhedron.icosahedron(fold = F)

# Defino los puntos
puntos = p.get_points()
edges_verts = fc.edges(puntos, M)
np.savetxt('edges.txt',edges_verts)

# Defino las caras
faces = fc.triangulos(p.get_points())

# Creo los puntos sobre las caras
puntos_face = []
for i in range(len(faces)):
    puntos_face.extend(fc.planos(faces[i],M))

np.savetxt('icosahedron.txt',puntos_face)
# Plot del icosahedro completo
# x_plot   = []
# y_plot   = []
# z_plot_v = []
# for i in range(len(puntos_face)):
#     x_plot.append(puntos_face[i][0])
#     y_plot.append(puntos_face[i][1])
#     z_plot_v.append(puntos_face[i][2])
# 
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(x_plot,y_plot,z_plot_v)
# plt.show()


# Calculo el área de integración
superficie_temp = p.intersect(hc)

# Calculo el área del virus que esta interaccionando
area = p.area_ico(hc,a)
#print('El área interactuante del icosaedro es: ',area)


# Plot de los vértices y las uniones 
#p.plot(show_labels = True)

# Selecciono los puntos que me interesan (los que están por debajo de la altura de corte)
puntos_temp = fc.CorteZ(puntos_face,hc)


# Itero sobre  hv, theta, hm y k para encontrar los valores que minimizan el volumen
# hv: altura a la que se encuentra el virus con respecto al plano
# theta: ratación alrededor del eje z que se le aplica al virus
# hm: altura máxima de la superficie de la molécula
# k: número de onda (2*pi/longitud de onda) de la superficie molecular
#for hv in np.linverts[union[0]][0] + i*deltaXspace(0.01*a,2*a,N):
for hv in [1]:

    #for t in np.linspace(0,2*np.pi/5,N):
    for t in [0]:

        # Traslado y roto el icosahedro según hv y t
        puntos_hv_t = fc.Transf(puntos_temp,[0,0,hv],[0,0,t]) 
        
        # Roto el área de integración
        superficie = fc.rot2d(superficie_temp,0)
        
        # Genero el mesh del virus
        Zv = fc.SupV(puntos_hv_t,R,M)
        
        # Plot de la sección del icosahedro elegida
        #xv = np.linspace(-R, R, M)
        #yv = np.linspace(-R, R, M)
        #x_plot   = []
        #y_plot   = []
        #z_plot_v = []
        #for i in range(len(xv)):
        #    for j in range(len(yv)):
        #        x_plot.append(xv[i])
        #        y_plot.append(yv[j])
        #        z_plot_v.append(Zv[j,i])
        #        #z_plot_m.append(Zm[j,i])
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
        #ax.scatter(x_plot,y_plot,z_plot_v)
        #plt.show()
        
        # Paso ahora a introducir los parámetros de la superficie
        #for hm in np.linspace(0,2*a,N):
        for hm in [0.2]: 
            #for k in np.linspace(0,10,N):
            for k in [1]:   

                # Genero el mesh de la molécula
                Zm = fc.SupM(fc.funcM,R,M,k,hm)

                # Compruebo que no haya intersección entre el virus y la molécula
                if fc.Intersec(Zm,Zv,M,R,superficie):
                    print('F')
                    break
                
                # Calculo los volúmenes de la molécula y el virus
                Vm = fc.Vmol(fc.funcM,hm,k,R,M,superficie)
                Vv = fc.Vvir(Zv,M,R,F,hv,hc,a)
                
                # Calculo el volumen libre entre ambos 
                Vl = Vv - Vm

                # Si este volumen es menor que el que teníamos previamente, los resultados se actualizan
                if Vl < V:
                    V = Vl
                    Resul = [V,hv,t,hm,k]

# Printeo los resultados
print(Resul)

#Guardo los puntos del virus y de la molécula en un formato cómodo para plotear
xv = np.linspace(-R, R, M)
yv = np.linspace(-R, R, M)
x_plot   = []
y_plot   = []
z_plot_v = []
z_plot_m = []
for i in range(len(xv)):
    for j in range(len(yv)):
        x_plot.append(xv[i])
        y_plot.append(yv[j])
        z_plot_v.append(Zv[j,i])
        z_plot_m.append(Zm[j,i])

x_plot_p = np.array(x_plot)
y_plot_p = np.array(y_plot)
z_plot_p_v = np.array(z_plot_v)
z_plot_p_m = np.array(z_plot_m)

np.savetxt('x.txt', xv, delimiter=',')
np.savetxt('y.txt', yv, delimiter=',')
np.savetxt('z_v.txt', Zv, delimiter=' ')
np.savetxt('z_m.txt', Zm, delimiter=' ')

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter(x_plot,y_plot,z_plot_v)
#ax.scatter(x_plot,y_plot,z_plot_m)
#plt.show()
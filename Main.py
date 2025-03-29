from shapely.geometry import Polygon, Point
import numpy as np
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import matplotlib.pyplot as plt
import math as math

for F in [2,3,5]:
    # Defino los parámetros que voy a usar 
    a = 1                                                                   # Longitud de la arista del icosaedro
    #F = 2                                                                   # Tipo de simetría (X-fold)
    R = a*np.sin(2*np.pi/5)                                                 # Distancia centro-vértice en el icosaedro
    M = 100                                                                 # Número de puntos creados sobre las caras
    N = 10
    hc_ = np.array([a*np.sqrt(1-1/(4*np.sin(np.pi/5)*np.sin(np.pi/5))),     # Altura de corte para los cálculos (calculada desde 
                a*np.sin(np.pi/3)*np.sin(np.pi-np.arccos(-np.sqrt(5)/3)),   # el punto más bajo del icosahedro)
                a*np.cos(np.pi/5)/2]) - 0.00001*a                                                       
    hc = hc_[0]*(F==5) + hc_[1]*(F==3) + hc_[2]*(F==2)
    # Inicializo el volumen (comienza en un valor muy elevado para asegurar que no hay problemas más adelante)
    V = 1e10*a

    # Creo una lista que guardará los valores de los parámetros en cada iteración
    V_list = []

    # Creo el icosaedro 
    p = Polyhedron.icosahedron(fold = F)

    # Defino los puntos
    puntos = p.get_points()
    edges_verts = fc.edges(puntos, M)
    np.savetxt('edges_' + str(F) + '.txt',edges_verts)

    # Defino las caras
    faces = fc.triangulos(p.get_points())

    # Creo los puntos sobre las caras
    puntos_face = []
    for i in range(len(faces)):
        puntos_face.extend(fc.planos(faces[i],M))

    np.savetxt('icosahedron_' + str(F) + '.txt',puntos_face)
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
    np.savetxt('interaction_'  + str(F) + '.txt',puntos_temp)

    # Itero sobre  hv, theta, hm y k para encontrar los valores que minimizan el volumen
    # hv: altura a la que se encuentra el virus con respecto al plano
    # theta: ratación alrededor del eje z que se le aplica al virus
    # hm: altura máxima de la superficie de la molécula
    # k: número de onda (2*pi/longitud de onda) de la superficie molecular
    for hv in np.linspace(0.5,3*a,N):
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
            for hm in [0,0.5,2]: 
                for k in [0,1,4]:
                    if (hm == 0 and k != 0) or (hm != 0 and k == 0):
                        continue
                    # Genero el mesh de la molécula
                    Zm = fc.SupM(fc.funcM,R,M,k,hm)

                    # Compruebo que no haya intersección entre el virus y la molécula
                    if fc.Intersec(Zm,Zv,M,R,superficie):
                        break
                    
                    # Calculo los volúmenes de la molécula y el virus
                    Vm = fc.Vmol(fc.funcM,hm,k,R,M,superficie)
                    Vv = fc.Vvir(Zv,M,R,F,hv,hc,a)

                    # Calculo el volumen libre entre ambos 
                    Vl = Vv - Vm
                    V_list.append([Vl,hv,t,hm,k])
                    # Si este volumen es menor que el que teníamos previamente, los resultados se actualizan
                    if Vl < V:
                        V = Vl
                        Resul = [V,hv,t,hm,k]

    # Guardo el resultado al final de la lista para fácil acceso
    V_list.append(Resul)
    np.savetxt('Parameters_hv_'+str(F)+'_fold.txt',np.array(V_list))
    # Printeo el resultado
    print('El volumen mínimo obtenido es V = ','{:03.2f}'.format(Resul[0]),'a^3 y los parámetros correspondientes son:')
    print('Altura del virus sobre el plano: h_v = ', '{:03.2f}'.format(Resul[1]),' a^3')
    print('Ángulo de rotación alrededor del eje z: theta = ', '{:03.2f}'.format(math.degrees(Resul[2])),'º')
    print('Altura de la corrugación sobre el plano: h_m = ', '{:03.2f}'.format(Resul[3]), ' a^3')
    if Resul[4] != 0:
        print('Longitud de onda de la superficie corrugada: lambda = ', '{:06.2f}'.format(2*np.pi/Resul[4]),' a')
    else:
        print('Número de onda de la superficie corrugada: k = ', '{:03.2f}'.format(Resul[4]),' a')



    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.scatter(x_plot,y_plot,z_plot_v)
    #ax.scatter(x_plot,y_plot,z_plot_m)
    #plt.show()
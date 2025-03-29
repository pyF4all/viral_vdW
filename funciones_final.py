from shapely.geometry import Polygon, Point
import numpy as np
from numpy import ravel, sin, cos
import pandas as pd
import scipy.interpolate
import matplotlib.pyplot as plt 
from polyhedron_final import Polyhedron



# Función que rota la superficie (dada por un polígono de la libreria shapely) de integración en un ángulo t
def rot2d(figura,t):
    
    # figura: Polígono que representa el perímetro de la superficie
    # t: ángulo rotado
    
    # Obtengo las coordenadas de los vértices
    xx,yy = figura.exterior.coords.xy

    # Roto los puntos y devuelvo el nuevo polígono rotado
    puntos_temp = np.array([[xx[i],yy[i]] for i in range(len(xx)-1)])
    R = np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
    puntos_temp2 = (R @ puntos_temp.T).T
    x = puntos_temp2[:,0]
    y = puntos_temp2[:,1]
    puntos = Polygon(list(zip(x,y)))

    return puntos

# Función que devuelve las coordenadas de los  3 vértices que forman cada cara en formato lista
def triangulos(verts):
    
    faces = [
         # Caras de la pirámide pentagonal inferior
         [verts[9], verts[4], verts[10]],
         [verts[9], verts[4], verts[8]],
         [verts[9], verts[8], verts[7]],
         [verts[9], verts[7], verts[11]],
         [verts[9], verts[11], verts[10]],

         # Caras intermedias
         [verts[7], verts[6], verts[11]],
         [verts[1], verts[6], verts[11]],
         [verts[1], verts[10],verts[11]],
         [verts[1], verts[5], verts[10]],
         [verts[4], verts[5], verts[10]],
         [verts[0], verts[4], verts[5]],
         [verts[0], verts[4], verts[8]],
         [verts[0], verts[8], verts[3]],
         [verts[3], verts[8], verts[7]],
         [verts[3], verts[6], verts[7]],

         # Caras de la pirámide pentagonal superior
         [verts[2], verts[3], verts[6]],
         [verts[2], verts[6], verts[1]],
         [verts[2], verts[1], verts[5]],
         [verts[2], verts[5], verts[0]],
         [verts[2], verts[0], verts[3]],
        ]     
    
    return faces


# Función que calcula los puntos de los planos del icosahedro y los pasa a formato (x,y,z)
def planos(points,M):

    # points:   lista de listas [[x,y,z],[...],[...]] con los xyz de los 3 puntos de la cara
    # M:        número de puntos (en una dirección) del meshgrid creado 
    
    # Obtengo los puntos de la cara
    p0, p1, p2 = points
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    # Calculo analíticamente el plano que forman los 3 puntos
    ux, uy, uz  = [x1-x0, y1-y0, z1-z0]
    vx, vy, vz  = [x2-x0, y2-y0, z2-z0]

    normal   = np.array([uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx])
    nullcomp = np.abs(normal) < 1e-6

    point  = np.array(p0)
    d = -point.dot(normal)

    # Creo los puntos sobre el plano en el que está contenida la cara
    if (np.sum(nullcomp) <= 1) and (not nullcomp[2]):
        x = np.linspace(min(x0,x1,x2),max(x0,x1,x2),M)
        y = np.linspace(min(y0,y1,y2),max(y0,y1,y2),M)
        xx, yy = np.meshgrid(x,y)
        zz = (-normal[0]*xx - normal[1]*yy - d)/normal[2]
    elif (np.sum(nullcomp) == 1) and (nullcomp[2]):
        x = np.linspace(min(x0,x1,x2),max(x0,x1,x2),M)
        z = np.linspace(min(z0,z1,z2),max(z0,z1,z2),M)
        xx, zz = np.meshgrid(x,z)
        yy = (-normal[0]*xx - d)/normal[1]
    else:
        if not nullcomp[0]:
            xx = -d/normal[0] * np.ones((M,M))
            y = np.linspace(min(y0,y1,y2),max(y0,y1,y2),M)
            z = np.linspace(min(z0,z1,z2),max(z0,z1,z2),M)
            yy, zz = np.meshgrid(y,z)
        if not nullcomp[1]:
            x = np.linspace(min(x0,x1,x2),max(x0,x1,x2),M)
            yy = -d/normal[1] * np.ones((M,M))
            z = np.linspace(min(z0,z1,z2),max(z0,z1,z2),M)
            xx, zz = np.meshgrid(x,z)
        if not nullcomp[2]:
            x = np.linspace(min(x0,x1,x2),max(x0,x1,x2),M)
            y = np.linspace(min(y0,y1,y2),max(y0,y1,y2),M)
            zz = -d/normal[2] * np.ones((M,M))
            xx, yy = np.meshgrid(x,y)
    
    # Creo la cara utilizando solo los puntos del plano contenidos en la proyección (en plano X-Y) de dicha cara
    triangulo = Polygon([Point(x0,y0),Point(x1,y1),Point(x2,y2)])
    plane = [[xx[0,i],yy[j,0],zz[j,i]] for i in range(M) for j in range(M) if Point(xx[0,i],yy[j,0]).within(triangulo)]

    if nullcomp[2]:
        triangulo = Polygon([Point(x0,z0),Point(x1,z1),Point(x2,z2)])
        plane = [[xx[0,i],yy[j,i],zz[j,0]] for i in range(M) for j in range(M) if Point(xx[0,i],zz[j,0]).within(triangulo)]

    return plane


# Función que corta el icosahedro a la altura (hc) que se vaya a considerar para la interacción
def CorteZ(puntos,hc):

    XYZ = []
    for i in range(len(puntos)):
        if puntos[i][2] <= hc:
            XYZ.append(puntos[i])

    return XYZ


# Función que convierte los puntos del icosahedro en un mesh alineado con el de la superficie
def SupV(puntos,R,M):

    # Cargo los valores de x,y,z 
    puntos = np.array(puntos)
    x = puntos[:,0]
    y = puntos[:,1]
    z = puntos[:,2]

    # Creo los vectores del mesh (cuadrado de tamaño justo para contener el icosahedro proyectado)
    xv = np.linspace(-R, R, M)
    yv = np.linspace(-R, R, M)

    # Creo el mesh
    X,Y = np.meshgrid(xv, yv)

    # Obtengo el mesh interpolando los valores iniciales
    Zv = scipy.interpolate.griddata((x,y),z,(X,Y),method='linear')

    # Limpio la interpolación poniendo los NaN a 0 
    for i in range(M):
        for j in range(M):
            if np.isnan(Zv[j,i]): 
                Zv[j,i] = 0

    return Zv

# Función matemática que define la superficie de la molécula en función de la posición y los parámetros hm y k
def funcM(x,y,hm,k):
    return hm *( np.sin(k*x + 3*np.pi/2)+np.sin(k*y + 3*np.pi/2) +2)/4


# Función que crea el mesh de la molécula
def SupM(funcM,R,M,k,hm):

    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)
    X,Y = np.meshgrid(x,y)
    
    Zm = funcM(X,Y,hm,k)
    return Zm


# Función que compruebasi hay intersección entre la molécula y el virus en algún punto
def Intersec(Zmol,Zv,M,R,Superficie):
    Zinter = (Zv - Zmol)
    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)
    
    for i in range(M):
        for j in range(M):
            if Point(x[i],y[j]).within(Superficie) and Zv[j,i] !=0 and Zinter[j,i] < 0: 
                return True
    
    return False

# Función que calcula el volumen bajo la molécula mediante sumas de Riemann 
def Vmol(funcM,hm,k,R,M,Superficie):

    # Mismos vectores de posición que los del mesh de la molécula
    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)

    # Tamaño de los cuadrados del mesh 
    deltaX = 2*R/(M-1)
    deltaY = 2*R/(M-1)
    base = deltaX*deltaY
    
    # Inicializo el volumen
    Vtotal = 0

    # Itero los puntos y hago las sumas de Riemann
    for i in range(M-1):
        for j in range(M-1):
            # Creo el prisma de Riemann correspondiente a la celda (i,j)
            polR = Polygon([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])

            # El valor de la función en cada prisma se considera como el valor de la función en el punto central
            altf = funcM(x[i]+deltaX/2,y[j]+deltaY/2,hm,k)   
            
            # Si el prisma esta contenido en el área de integración se suma el volumen completo, si solo interseca 
            # se suma la mitad (para paliar el volumen que se desprecia debido a la forma de la superficie)
            check = polR.within(Superficie) + 0.5*(not polR.within(Superficie))*polR.intersects(Superficie)
            Vtotal += base*altf*check

    return Vtotal


# Función que calcula el volumen bajo el virus mediante sumas de Riemann
def Vvir(Zv,M,R,F,hv,hc,a):

    # En el caso 5 fold podemos hacerlo geométricamente, ahorrando tiempo
    if F == 5:
        # Calculamos los parámetros necesarios
        L = a
        ap = a/(2*np.tan(np.pi/5))
    
        # Volumen del prisma menos volumen de la pirámide
        V = 5*L*ap*(hc/3+hv/2)
        return V

    # Tamaño de los cuadrados del mesh
    deltaX = 2*R/(M-1)
    deltaY = 2*R/(M-1)
    base = deltaX*deltaY

    # Inicializo el volumen
    Vtotal = 0

    # Itero sobre todo el mesh haciendo las sumas de Riemann
    for i in range(M-1):
        for j in range(M-1):

            # El valor de la función en cada prisma se considera como el valor máximo de la función en los 4 puntos
            altf = max([Zv[j,i], Zv[j+1,i], Zv[j,i+1], Zv[j+1,i+1]])
            Vtotal += base*altf

    return Vtotal

# Create edges
def edges(verts,M):
    
    edges = []
    # Defino las aristas
    unions = [
        (1,9),(1,3),(1,4),(1,5),(1,6),
        (2,3),(2,6),(2,7),(2,11),
        (3,7),(3,6),(3,4),
        (4,8),(4,9),
        (5,11),(5,9),(5,10),(5,6),
        (6,11),
        (4,7),(7,12),(7,8),
        (8,12),(8,9),(8,10),
        (9,10),
        (10,11),(10,12),
        (11,12),
        (2,12)
    ]    

    for union in unions:
        deltaX = (verts[union[1]-1][0] - verts[union[0]-1][0])/M
        deltaY = (verts[union[1]-1][1] - verts[union[0]-1][1])/M
        deltaZ = (verts[union[1]-1][2] - verts[union[0]-1][2])/M
        
        for i in range(M):
            edges.append([verts[union[0]-1][0] + i*deltaX, verts[union[0]-1][1] + i*deltaY, verts[union[0]-1][2] + i*deltaZ])

    return edges


# Function that calculates the average distance between the particles forming the virus and a plane situated at a distance D
def dist_med(F, a, D, k, hm, hc):
    
    # F:  fold
    # D:  minimum distance between icosahedron and plane
    # a:  length of icosahedron edge
    # k:  wave number of surface
    # hm: max height of surface
    # hc: maximum height considered (if it's more that the one allowed by the script it becomes that, the allowed heights are
    # described in hc_max)


    R = a*np.sin(2*np.pi/5)
    M = 100
    hc_max =np.empty(3)
    hc_max[0] = a*np.sqrt(1-1/(4*np.sin(np.pi/5)*np.sin(np.pi/5)))          # 5-fold
    hc_max[1] = a*np.sin(np.pi/3)*np.sin(np.pi-np.arccos(-np.sqrt(5)/3))    # 3-fold
    hc_max[2] =  a*np.cos(np.pi/5)/2                                        # 2-fold
    
    hc_max += -0.001 # To make sure we don't run into computational errors
    selector = 0*(F==5) + 1*(F==3) + 2*(F==2) # Selects the corresponding fold

    hc_used = min(hc,hc_max[selector])

    # Creo el icosaedro 
    p = Polyhedron.icosahedron(fold = F)

    # Defino los puntos
    puntos = p.get_points()
    
    # Defino las caras
    faces = triangulos(p.get_points())

    # Creo los puntos sobre las caras
    puntos_face = []
    for i in range(len(faces)):
        puntos_face.extend(planos(faces[i],M))
    
    # Selecciono los puntos que me interesan (los que están por debajo de la altura de corte)
    puntos_temp = CorteZ(puntos_face,hc_used)

    puntos = Transf(puntos_temp,[0,0,D],[0,0,0])
    supV = SupV(puntos,R,M)
    supM = SupM(funcM,R,M,k,hm)

    distances = supV - supM 
    d_m = np.average(distances)

    return d_m











######## FUNCIONES EN DESUSO POR AHORA ########

# Función que crea un pentágono regular con radio centro-vértice R y rotado un ángulo t con respecto al eje Y
def pentagono(R,t):
    
    c1 = R*np.cos(2*np.pi/5)
    c2 = R*np.cos(np.pi/5)
    s1 = R*np.sin(2*np.pi/5)
    s2 = R*np.sin(4*np.pi/5)

    P0 = Point(-R*np.sin(t),R*np.cos(t))
    P1 = Point(s1*np.cos(t)-c1*np.sin(t),s1*np.sin(t)+c1*np.cos(t))
    P2 = Point(s2*np.cos(t)+c2*np.sin(t),s2*np.sin(t)-c2*np.cos(t))
    P3 = Point(-s2*np.cos(t)+c2*np.sin(t),-s2*np.sin(t)-c2*np.cos(t))
    P4 = Point(-s1*np.cos(t)-c1*np.sin(t),-s1*np.sin(t)+c1*np.cos(t)) 

    return Polygon([P0,P1,P2,P3,P4])

# Función que crea los vértices del icosahedro
def vertices(a):

    def vertex(x, y, z):

        return [ i*a/2  for i in (x,y,z)]
    
    PHI = (1 + np.sqrt(5)) / 2

    verts = [

            vertex(-1,  PHI, 0),
            vertex( 1,  PHI, 0),
            vertex(-1, -PHI, 0),
            vertex( 1, -PHI, 0),

            vertex(0, -1, PHI),
            vertex(0,  1, PHI),
            vertex(0, -1, -PHI),
            vertex(0,  1, -PHI),

            vertex( PHI, 0, -1),
            vertex( PHI, 0,  1),
            vertex(-PHI, 0, -1),
            vertex(-PHI, 0,  1),            
        
            ]  
    return verts


# Función que devuelve los índices de los puntos que forman cada cara
def faces():
        faces = [
         # 5 faces around point 0
         (1,12,6),
         (1,6,2),
         (1,2,8),
         (1,8,11),
         (1,11,12),

         # Adjacent faces
         (2,6,10),
         (6,12,5),
         (12,11,3),
         (11,8,7),
         (8,2,9),

         # 5 faces around 3
         (4,10,5),
         (4,5,3),
         (4,3,7),
         (4,7,9),
         (4,9,10),

         # Adjacent faces
         (5,10,6),
         (3,5,12),
         (7,3,11),
         (9,7,8),
         (10,9,2)
        ]   

        return faces

# Función que da las rotaciones y traslaciones correspondientes en función de la simetría 
# S = 0 da 5fold, S = 1 da 3fold, S = 2 da 2fold
def Symmetry(S,a):

    fold_5 = [np.pi/6,0,0]
    fold_3 = []
    fold_2 = [0,0,0] 

    traslacion = [
        a*np.sin(2*np.pi/5),
        np.sqrt(3)/12*(3+np.sqrt(5)),
        a*np.cos(np.pi/5)]        

    return (S == 0)*fold_5 +(S == 1)*fold_3 + (S == 2)*fold_2, [0,0,traslacion[S]]

# Función que aplica traslaciones y/o rotaciones (en sentido antihorario) a los puntos del icosahedro
def Transf(puntos,Traslacion,Rotaciones):

    # Traslación es una lista con los valores de [x,y,z] a trasladar
    # Rotaciones es una lista con los ángulos a rotar en cada eje [theta(X),theta(Y),theta(Z)]
    
    # Hago una traslación a todos los puntos que correspondan
    Traslacion = np.array(Traslacion)
    puntos = np.array(puntos)

    

    # Defino la matriz de rotación   
    a = Rotaciones[0]
    b = Rotaciones[1]
    c = Rotaciones[2]
    R =np.array(
        [[cos(b)*cos(c), sin(a)*sin(b)*cos(c)-cos(a)*sin(c), cos(a)*sin(b)*cos(c)+sin(a)*sin(c)],
        [cos(b)*sin(c), sin(a)*sin(b)*sin(c)+cos(a)*cos(c), cos(a)*sin(b)*sin(c)-sin(a)*cos(c)],
        [-sin(b), sin(a)*cos(b), cos(a)*cos(b)]
        ])
    
    # Aplico la rotación sobre puntos
    puntos = (R @ puntos.T).T + Traslacion

    return puntos 


def Vvir_antiguo(Zv,N,R,S,Superficie,hv,hc,a):

    # En el caso 5 fold podemos hacerlo geométricamente, ahorrando tiempo
    if S == 5:
        # Calculamos los parámetros necesarios
        L = a
        ap = a/(2*np.tan(np.pi/5))
    
        # Volumen del prisma menos volumen de la pirámide
        V = 5*L*ap*(hc/3+hv/2)
        return V


    # Mismos vectores de posición que los del mesh de la molécula
    x = np.linspace(-R,R,N)
    y = np.linspace(-R,R,N)

    # Tamaño de los cuadrados del mesh
    deltaX = 2*R/(N-1)
    deltaY = 2*R/(N-1)
    base = deltaX*deltaY

    # Inicializo el volumen
    Vtotal = 0

    # Itero sobre todo el mesh haciendo las sumas de Riemann
    for i in range(N-1):
        for j in range(N-1):
            # Creo la base del prisma de Riemann correspondiente a la celda (i,j)
            polR = Polygon([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])

            # El valor de la función en cada prisma se considera como el valor medio de la función en los 4 puntos
            altf = max([Zv[j,i], Zv[j+1,i], Zv[j,i+1], Zv[j+1,i+1]])
            
            # Si el prisma esta contenido en el área de integración se suma el volumen completo, si solo interseca 
            # se suma la mitad (para paliar el volumen que se desprecia debido a la forma de la superficie)
            check = polR.within(Superficie) + 0.5*(not polR.within(Superficie))*polR.intersects(Superficie)
            Vtotal += base*altf#*check

    return Vtotal

# Función que calcula el volumen libre entre el prisma pentagonal (simetría 5-fold) y la superficie 
def Vol5fold(a,hv,Vmol):

    # Calculamos los parámetros necesarios
    L = a
    ap = a/(2*np.tan(np.pi/5))
    h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))

    # Volumen del prisma menos volumen de la pirámide
    V = 5*L*ap*(h/3+hv/2)

    # Obtenemos el volumen libre
    Vlibre = V - Vmol

    return Vlibre    
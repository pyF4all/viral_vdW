from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random as random




# We need a function that calculates the derjaguin approximation to the vdW force for two equal paralel polygons
def derjaguin_polygon(area, D, d, A_H):

    # area: area of the polygon
    # D: distance between both polygons
    # d: depth considered for each polygon faced plane
    # A_H: Hammaker constant

    # Utilizamos directamente la solución analítica ya que no va a variar entre casos distintos, y solo se ve modificado 
    # por un factor correspondiente al area

    F = area*A_H/(6*np.pi)*(2/(D+d)**3-1/D**3-1/(D+2*d)**3)
    return F


# We now create a function that projects our force onto the (x,y,z) coordinates 
def projection(normal):
    # normal: vector perpendicular (and outward facing) to the surface

    ux = np.array([1,0,0])
    uy = np.array([0,1,0])
    uz = np.array([0,0,1])

    x = np.dot(normal,ux)*ux
    y = np.dot(normal,uy)*uy
    z = np.dot(normal,uz)*uz

    ur = np.array([x, y, z])

    return ur

# Function
def vdW_face(fold, a, hc, D, d, A_H):

    ico = Polyhedron.icosahedron(fold=fold)
    areas, normals = Polyhedron.area_normal(ico, hc, a)

    areas = np.array(areas)
    normals = np.array(normals)
    F_temp = np.empty([len(areas),3])

    for i in range(len(areas)):
        F_module = abs(derjaguin_polygon(areas[i], D, d, A_H))
        #print(projection(normals[i])[2,:])
        proj = np.array(projection(normals[i])[2,:])
        F_temp[i,:] = np.array(F_module*proj)

    F = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
    return F

def F_5_fold(D,d,A_H,N):

    area = 5*np.sqrt(3)/4
    F_total = derjaguin_polygon(area, D, d, A_H)
    return F_total


# The legend for the d values must be changed manually
def plot_fold(fold, a, D, d, hc, A_H):

    F = np.empty([len(d), len(D)])

    for i in range(len(d)):
        for j in range(len(D)):
            F[i,j] = vdW_face(fold, a, hc, D[j], d[i], A_H)[2]

    for i in range(len(d)):
        plt.plot(D, F[i,:], '.')
    if fold ==5:
        plt.title('FvdW 5fold-5fold for different plane depth')
    elif fold ==3:
        plt.title('FvdW 3fold-3fold for different plane depth')
    elif fold ==2:
        plt.title('FvdW 2fold-2fold for different plane depth')

    plt.xlabel('Distancia mínima (a)')
    plt.ylabel('Fuerza de van der Waals ($A_H$/a)')
    plt.title('Fuerza vdW frente a profundidad de plano ($h_c =0.5 a$)')
    plt.grid()
    plt.legend(['d=0.01','d=0.1','d=0.2','d=0.3','d=0.4','d=0.5','d=10'])
    plt.savefig('d.png')


def plot_folds(a, D, d, hc, A_H):
    
    F5 = np.empty(len(D))
    F3 = np.empty(len(D))
    F2 = np.empty(len(D))

    for j in range(len(D)):

        F5[j] = vdW_face(5, a, hc, D[j], d, A_H)[2]
        F3[j] = vdW_face(3, a, hc, D[j], d, A_H)[2]
        F2[j] = vdW_face(2, a, hc, D[j], d, A_H)[2]

    plt.plot(D, F5,'.')
    plt.plot(D, F3,'.')
    plt.plot(D, F2,'.')
    plt.title('Fvdw icosahedral surfaces for each fold')
    plt.legend(['5-fold','3-fold','2-fold'])
    plt.xlabel('Minimum distance (D)')
    plt.ylabel('Van der Waals force in Z axis')
    plt.show()

def plot_comp(a, D, d, A_H):
    
    F5 = np.empty(len(D))
    F3 = np.empty(len(D))
    F2 = np.empty(len(D))

    hc5 = 2*a*np.sin(2*np.pi/5)
    hc3 = 2*np.sqrt(3)/12*(3+np.sqrt(5))*a
    hc2 = 2*a*np.cos(np.pi/5)

    for j in range(len(D)):

        F5[j] = vdW_face(5, a, hc5, D[j], d, A_H)[2]
        F3[j] = vdW_face(3, a, hc3, D[j], d, A_H)[2]
        F2[j] = vdW_face(2, a, hc2, D[j], d, A_H)[2]

    plt.plot(D, F5,'.')
    plt.plot(D, F3,'.')
    plt.plot(D, F2,'.')
    plt.title('Fvdw icosahedral surfaces for each fold (complete)')
    plt.legend(['5-fold','3-fold','2-fold'])
    plt.xlabel('Minimum distance (D')
    plt.ylabel('Van der Waals force in Z axis')
    plt.show()

def plot_hc(a, D, d, A_H, N, axis):
    
    F5 = np.empty(N)
    F3 = np.empty(N)
    F2 = np.empty(N)

    hc5 = np.linspace(0.001, 2*a*np.sin(2*np.pi/5), N)
    hc3 = np.linspace(0.001, 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, N)
    hc2 = np.linspace(0.001, 2*a*np.cos(np.pi/5), N)

    for j in range(len(hc5)):

        F5[j] = vdW_face(5, a, hc5[j], D, d, A_H)[axis]
        F3[j] = vdW_face(3, a, hc3[j], D, d, A_H)[axis]
        F2[j] = vdW_face(2, a, hc2[j], D, d, A_H)[axis]

    F5_plot = F5[np.argmin(F5)]
    F3_plot = F3[np.argmin(F3)]
    F2_plot = F2[np.argmin(F2)]

    plt.plot(hc2, F2,'.')
    plt.plot(hc3, F3,'.')
    plt.plot(hc5, F5,'.')
    plt.plot(hc2[np.argmin(F2)],F2_plot,'r.',ms=12)
    plt.plot(hc3[np.argmin(F3)],F3_plot,'r.',ms=12)
    plt.plot(hc5[np.argmin(F5)],F5_plot,'r.',ms=12)
    plt.title('Fvdw frente a altura considerada en superficie corrugada')
    plt.legend(['2-fold','3-fold','5-fold','Máximos'])
    plt.xlabel('Altura considerada (a)')
    if axis == 0:
        plt.ylabel('Van der Waals force in X axis')
    if axis == 1:
        plt.ylabel('Van der Waals force in Y axis')
    if axis == 2:
        plt.ylabel('Fuerza de van der Waals ($\mathrm{A_H}$/a)')    
    plt.grid()
    plt.show()
    #plt.savefig('corrug_hc.png')

# Function to plot vectors in 3D with plotly
def vector_plot(tvects,is_vect=True,orig=[0,0,0]):
    """Plot vectors using plotly"""

    if is_vect:
        if not hasattr(orig[0],"__iter__"):
            coords = [[orig,np.sum([orig,v],axis=0)] for v in tvects]
        else:
            coords = [[o,np.sum([o,v],axis=0)] for o,v in zip(orig,tvects)]
    else:
        coords = tvects

    data = []
    for i,c in enumerate(coords):
        X1, Y1, Z1 = zip(c[0])
        X2, Y2, Z2 = zip(c[1])
        vector = go.Scatter3d(x = [X1[0],X2[0]],
                              y = [Y1[0],Y2[0]],
                              z = [Z1[0],Z2[0]],
                              marker = dict(size = [0,5],
                                            color = ['blue'],
                                            line=dict(width=20,
                                                      color='DarkSlateGrey')),
                              name = 'hc'+str(i+1))
        data.append(vector)

    layout = go.Layout(
             margin = dict(l = 4,
                           r = 4,
                           b = 4,
                           t = 4)
                  )
    fig = go.Figure(data=data,layout=layout)
    fig.show()


# Force vector vs hc
def plot_force(a, A_H, D, d, N, fold):

    F = []
    select = int(0 + 1*(fold==3) + 2*fold==(3))
    hc_temp = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)]) - 0.001
    hc = np.linspace(0.001, hc_temp[select], N)

    for j in range(len(hc)):

        F_temp = vdW_face(fold, a, hc[j], D, d, A_H)
        F.append(F_temp/np.linalg.norm(F_temp))

    vector_plot(F)



# We are now going to add a series of functions that will allow for the simulation of simplified "docking proteins" that 
# modify the vdw foce

# Function that calculates the vdw force between a semi-sphere and a plane. We are using the analytical result directly.
def vdw_semi(r, A_H, D_r):

    F = -A_H/6 * (1/(r+D_r) + (r-D_r)/(2*D_r**2))

    return F

# Now we modify the vdw_face function to incorporate this new addition. This calculation only works for faces with full 
# (5-fold y solo para piramide pentagonal por ahora)
def vdW_face_r(fold, a, r, hc, D, d, A_H):

    D_r = D-r
    ico = Polyhedron.icosahedron(fold=fold)
    areas, normals = Polyhedron.area_normal(ico, hc, a)

    areas = np.array(areas)
    normals = np.array(normals)
    F_temp = np.empty([len(areas),3])
    F_temp_r = np.empty([len(areas), 3])

    F_module_r = abs(vdw_semi(r,A_H,D_r))

    for i in range(len(areas)):
        F_module = abs(derjaguin_polygon(areas[i], D, d, A_H))
        #print(projection(normals[i])[2,:])
        proj = np.array(projection(normals[i])[2,:])
        F_temp[i,:] = np.array(F_module*proj)
        F_temp_r[i,:] = np.array(F_module_r*proj)

    F = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
    F_r = np.array([np.sum(F_temp_r[:,0]), np.sum(F_temp_r[:,1]), np.sum(F_temp_r[:,2])])    
    
    Ftotal = F+F_r

    return Ftotal


def vdW_face_r_good(fold, a, r, hc, D, d, A_H):

    D_r = D-r
    ico = Polyhedron.icosahedron(fold=fold)
    areas, normals = Polyhedron.area_normal(ico, hc, a)
    centers = Polyhedron.face_center(ico)
    
    areas = np.array(areas)
    normals = np.array(normals)
    
    F_temp = np.empty([len(areas),3])
    F_temp_r = np.empty([len(areas), 3])

    F_module_r = abs(vdw_semi(r,A_H,D_r))

    for i in range(len(areas)):
        F_module = abs(derjaguin_polygon(areas[i], D, d, A_H))
        #print(projection(normals[i])[2,:])
        proj = np.array(projection(normals[i])[2,:])
        F_temp[i,:] = np.array(F_module*proj)

        if (hc > (centers[i][2] + r/2)): 
            F_temp_r[i,:] = np.array(F_module_r*proj)
        else:
            F_temp_r[i,:] = np.array(0*proj)
    F = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
    F_r = np.array([np.sum(F_temp_r[:,0]), np.sum(F_temp_r[:,1]), np.sum(F_temp_r[:,2])])
    #print(F_r)    
    
    Ftotal = F+F_r

    return Ftotal

# Returns the relation F_mod/F
def difference_protein(fold, a, r, hc, D, d, A_H):

    F     = vdW_face(fold,a,hc,D,d,A_H)[2]
    F_mod = vdW_face_r_good(fold,a,r,hc,D,d,A_H)[2]

    difference = abs(F_mod/F)

    return difference

def plot_folds(a, D, d, hc, A_H):
    
    F5 = np.empty(len(D))
    F3 = np.empty(len(D))
    F2 = np.empty(len(D))

    for j in range(len(D)):

        F5[j] = vdW_face(5, a, hc, D[j], d, A_H)[2]
        F3[j] = vdW_face(3, a, hc, D[j], d, A_H)[2]
        F2[j] = vdW_face(2, a, hc, D[j], d, A_H)[2]

    plt.plot(D, F2,'.')
    plt.plot(D, F3,'.')
    plt.plot(D, F5,'.')
    plt.title('Fuerza de van der Waals entre cada simetría y la superficie corrugada')
    plt.legend(['2-fold','3-fold','5-fold'])
    plt.xlabel('Distancia mínima (a)')
    plt.ylabel('Fuerza de interacción vdW ($\mathrm{A_H}$/a)')
    plt.yscale("symlog")
    plt.grid()
    plt.show()
    #plt.savefig('corrug_folds.png')

def plot_protein_fold(fold,a,hc,D,d,A_H):

    r = [0.01, 0.02, 0.045]
    F = np.empty(len(D))
    F_mod = np.empty(len(D))
    F_mod2 = np.empty(len(D))
    F_mod3 = np.empty(len(D))

    for i in range(len(D)):

        F[i] = vdW_face(fold,a,hc,D[i],d,A_H)[2]
        F_mod[i] = vdW_face_r_good(fold,a,r[0],hc,D[i],d,A_H)[2]
        F_mod2[i] = vdW_face_r_good(fold,a,r[1],hc,D[i],d,A_H)[2]
        F_mod3[i] = vdW_face_r_good(fold,a,r[2],hc,D[i],d,A_H)[2]

    plt.plot(D,F,'.')
    plt.plot(D,F_mod,'.')
    plt.plot(D,F_mod2,'.')
    plt.plot(D,F_mod3,'.')
    plt.title('Fvdw 5-fold with proteins for hc = 0.5')
    plt.legend(['5-fold','r = 0.01','r=0.02','r=0.045'])
    plt.xlabel('Minimum distance (D)')
    plt.ylabel('Van der Waals force in Z axis')
    plt.grid()
    plt.show()



def plot_protein_hc(fold,a,D,d,r,A_H):

    hc = [0.1, 0.5, 1, 1.5]
    F = np.empty(len(D))
    F2 = np.empty(len(D))
    F3 = np.empty(len(D))
    F4 = np.empty(len(D))
    F_mod = np.empty(len(D))
    F_mod2 = np.empty(len(D))
    F_mod3 = np.empty(len(D))
    F_mod4 = np.empty(len(D))

    for i in range(len(D)):

        F[i] = vdW_face(fold,a,hc[0],D[i],d,A_H)[2]
        F_mod[i] = vdW_face_r_good(fold,a,r,hc[0],D[i],d,A_H)[2]
        
        F2[i] = vdW_face(fold,a,hc[1],D[i],d,A_H)[2]
        F_mod2[i] = vdW_face_r_good(fold,a,r,hc[1],D[i],d,A_H)[2]

        F3[i] = vdW_face(fold,a,hc[2],D[i],d,A_H)[2]
        F_mod3[i] = vdW_face_r_good(fold,a,r,hc[2],D[i],d,A_H)[2]

        F4[i] = vdW_face(fold,a,hc[3],D[i],d,A_H)[2]
        F_mod4[i] = vdW_face_r_good(fold,a,r,hc[3],D[i],d,A_H)[2]

    plt.plot(D,F,'.')
    plt.plot(D,F_mod,'.')
    plt.plot(D,F2,'.')
    plt.plot(D,F_mod2,'.')
    plt.plot(D,F3,'.')
    plt.plot(D,F_mod3,'.')
    plt.plot(D,F4,'.')
    plt.plot(D,F_mod4,'.')

    plt.title('Fvdw 5-fold with proteins for hc = 0.5')
    plt.legend(['hc = 0.1','hc = 0.1 (mod)','hc = 0.5','hc = 0.5 (mod)','hc = 1','hc = 1 (mod)','hc = 1.5','hc = 1.5 (mod)'])
    plt.xlabel('Minimum distance (D)')
    plt.ylabel('Van der Waals force in Z axis')
    plt.grid()
    plt.show()
    

def difference_r(fold, a, r, D, d, A_H):

    hc = [0.1, 0.5, 1, 1.5]
    dif  = np.empty(len(r))
    dif2 = np.empty(len(r))
    dif3 = np.empty(len(r))
    dif4 = np.empty(len(r))

    for i in range(len(r)):
        dif[i]  = difference_protein(fold,a,r[i],hc[0],D,d,A_H)
        dif2[i] = difference_protein(fold,a,r[i],hc[1],D,d,A_H)
        dif3[i] = difference_protein(fold,a,r[i],hc[2],D,d,A_H)
        dif4[i] = difference_protein(fold,a,r[i],hc[3],D,d,A_H)
    
    plt.plot(r,dif,'.')
    plt.plot(r,dif2,'.')
    plt.plot(r,dif3,'.')
    plt.plot(r,dif4,'.')
    plt.xlabel('Proteic radious (r)')
    plt.ylabel('F_protein/F')
    plt.title('Ratio of vdW forces with/without proteins (D = 0.06)')
    plt.grid()
    plt.legend(['hc = 0.1', 'hc = 0.5', 'hc = 1', 'hc = 1.5'])
    plt.show()

#r =0.03
#plot_protein_hc(5,a,D,d,r,A_H)
#a = 1
#A_H = 1
#D = 0.06
#d = 0.1
#hc = 0.5
#N = 5
#fold = 3
#r = np.linspace(0.01,0.055)
#difference_r(fold,a,r,D,d,A_H)
#ico = Polyhedron.icosahedron(fold=5, rot_extra=[45*180/np.pi,0,0])
#Polyhedron.plot(ico)




### FUNCIONES PARA  CALCULAR LAS INTERACCIONES DE VDW CON PROTEÍNAS EN DETERMINADAS POSICIONES ###
##### POSIBLE MEJORA: CALCULAR UNA VEZ LOS VECTORES Y LAS Z PARA LA POSICIÓN BASE EN LA QUE SE CREA EL ICOSAEDRO, 
##### Y QUE LA FUNCIÓN LO ÚNICO QUE HAGA SEA ROTAR, TRASLADAR Y SELECCIONAR LOS PUNTOS
 
# Function that calculates and normalizes the vector from the center of the icosahedron to a certain point
def vector(center, point):

    out_temp = np.array([point[0]-center[0], point[1]-center[1], point[2]-center[2]])
    out_norm = out_temp/np.linalg.norm(out_temp)

    return out_norm

# function that calculates the middle point of two other 3d points
def middle(p1,p2):

    out = [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2]

    return out


# Function that calculates the closest point between a given reference and a list of points and returns its index
def closest(reference, points):
    # reference: point from which we are calculating the distance [x,y,z]
    # list: list of points we are checking [[x,y,z], [x,y,z], ...]
    
    reference = np.array(reference)
    points = np.array(points)

    distances = np.sum((points-reference)**2)
    index = np.argmin(distances)

    return index


# Función que calcula las posibles posiciones de las proteínas sobre cada cara, returning the vectors conecting them
# with the icosahedron's center, and their z coordinate
def proteins(ico, number):

    # ico:    icosahedron from the polyhedron class
    # number: number of positions to be considered (1, 3, 4, 7, 10 or 13)

    # List containing the [x,y,x] coordinates of the icosahedron's center
    ico_center = Polyhedron.center(ico)

    # List containing the coordinates of all vertex, grouped on sets of 3 corresponding to each face
    vertex  = Polyhedron.face_vertices(ico)

    # List containing the coordinates of the center of each face
    centers = Polyhedron.face_center(ico)

    # List that will contain the results, the first dimension corresponds with the face in which we are, and the second one 
    # contains the vectors and the z coordinate of their points (list size is 20xnumber)
    out = [[0] * number for i in range(20)]

    # Calculate the first set of points and their normals (the centers of each face)
    for i in range(len(centers)):
            normal = vector(ico_center, centers[i])
            out[i][0] = [normal,centers[i][2]]
    

    # If we only want the centers, we can return out now
    if number == 1: 
        return out
    
    # Now we calculate the positions of the points that lie in the middle of the center of each face and the vertices
    # and the corresponding normals
    p3_1 = []
    p3_2 = []
    p3_3 = []
    for i in range(len(vertex)):
        p3_1.append([middle(vertex[i][0], centers[i])])
        p3_2.append([middle(vertex[i][1], centers[i])])
        p3_3.append([middle(vertex[i][2], centers[i])])

        normal3_1 = vector(ico_center, p3_1[i])
        normal3_2 = vector(ico_center, p3_2[i])
        normal3_3 = vector(ico_center, p3_3[i])
        
        out[i][1] = [normal3_1, p3_1[i][2]]
        out[i][2] = [normal3_2, p3_2[i][2]]
        out[i][3] = [normal3_3, p3_3[i][2]]

    
    # Return 4 points
    if number == 4:
        return out
    
    # Next we have to calculate the points situated in the middle of the p3 points
    p7_1 = []
    p7_2 = []
    p7_3 = []

    for i in range(len(vertex)):

        p7_1.append([middle(p3_1[i],p3_2[i])])
        p7_2.append([middle(p3_2[i],p3_3[i])])
        p7_3.append([middle(p3_3[i],p3_1[i])])
        
        normal7_1 = vector(ico_center, p7_1[i])
        normal7_2 = vector(ico_center, p7_2[i])
        normal7_3 = vector(ico_center, p7_3[i])

        out[i][4] = [normal7_1, p7_1[i][2]]
        out[i][5] = [normal7_2, p7_2[i][2]]
        out[i][6] = [normal7_3, p7_3[i][2]]

    # If we want 7 points
    if number == 7:
        return out

    # Points in the middle of p3 and the vertex
    p10_1 = []
    p10_2 = []
    p10_3 = []

    for i in range(len(vertex)):

        p10_1.append([middle(p3_1[i],vertex[i][closest(p3_1[i],[vertex[i][0],vertex[i][0],vertex[i][2]])])])
        p10_2.append([middle(p3_2[i],vertex[i][closest(p3_2[i],[vertex[i][0],vertex[i][0],vertex[i][2]])])])
        p10_3.append([middle(p3_3[i],vertex[i][closest(p3_3[i],[vertex[i][0],vertex[i][0],vertex[i][2]])])])
        
        normal10_1 = vector(ico_center, p10_1[i])
        normal10_2 = vector(ico_center, p10_2[i])
        normal10_3 = vector(ico_center, p10_3[i])

        out[i][7] = [normal10_1, p10_1[i][2]]
        out[i][8] = [normal10_2, p10_2[i][2]]
        out[i][9] = [normal10_3, p10_3[i][2]]


    # If we want 10 points
    if number == 10:
        return out


# Finally, we add the points situated in the middle of each pair of p10_X
    p13_1 = []
    p13_2 = []
    p13_3 = []

    for i in range(len(vertex)):

        p13_1.append([middle(p10_1[i],p10_2[i])])
        p13_2.append([middle(p10_2[i],p10_3[i])])
        p13_3.append([middle(p10_3[i],p10_1[i])])
        
        normal13_1 = vector(ico_center, p13_1[i])
        normal13_2 = vector(ico_center, p13_2[i])
        normal13_3 = vector(ico_center, p13_3[i])

        out[i][10] = [normal13_1, p13_1[i][2]]
        out[i][11] = [normal13_2, p13_2[i][2]]
        out[i][12] = [normal13_3, p13_3[i][2]]

    # If we want 13 points 
    return out


# Now that we can calculate the possible positions for points in each face, we need a funcion that calculates the vdw force
# exherted by semispheres of radious r situated in said positions. We will create two functions, one where all faces have 
# the same number of "ocupied points" and one where said number is randomised

def symmetric_complete_vdw(ico, r, D, A_H, N, hc):

    # ico: icosahedron
    # N:   total number of points per face
    
    
    # First we calculate the possible points and their normal vectors
    points = proteins(ico,N)
    
    # Now we calculate the vdW force corresponding to a semi sphere of the chosen radious
    F_module = vdw_semi(r,A_H,D-r)

    # Now we iterate over all points and add the force corresponding to the ones with a z coordinate smaller than hc
    F = np.array([0,0,0])
    for i in range(20):
        for j in range(N):
            F += F_module * points[i][j][0] * ((points[i][j][1] + r/2) > hc)
    
    return F

# Function that calculates the vdW force generated in the interaction between a number n of semi spheres situated in each
# face of an icosahedron when interacting with a planar surface opposite to said face.
def point_vdw(ico, r, D, A_H, N, n, hc, sym):

    # ico: icosahedron
    # N:   total number of points per face
    # n:   number of points to be "filled" with half-spheres in each face (can be a list if it varies from face to face in
    #      the assymetrical case)
    # sym: True if we want all faces to have the same distribution of points in the faces, False if we want the distribution
    #      to be random and vary from face to face
    
    # First we consider if n is a list or an integer
    if type(n) is not list:
        n = n*np.ones(20, dtype=int)

    # After that we calculate the possible points and their normal vectors
    points = proteins(ico,N)
    
    # Now we calculate the vdW force corresponding to a semi sphere of the chosen radious
    F_module = vdw_semi(r,A_H,D-r)

    # Define the force vector
    F = np.array([0,0,0])

    # Symmetrical case
    if sym == True:

        # Now we generate the points that will be taken into consideration
        points_considered =random.sample(range(N), n[0])

        # We iterate over all considered points and add the force corresponding to the ones with a z coordinate smaller 
        # than hc
        for i in range(20):
            for j in points_considered:
                F = F + F_module * points[i][j][0] * ((points[i][j][1] + r/2) < hc)

        return F
    
    # Assymetrical case
    elif sym == False:

        # Now we iterate over all considered points and add the force corresponding to the ones with a z coordinate smaller 
        # than hc
        for i in range(20):
            # For each face generate a different set of points to take into consideration
            points_considered =random.sample(range(0,N), n[i])
            for j in points_considered:
                F = F + F_module * points[i][j][0] * ((points[i][j][1] + r/2) < hc)

        return F



D = np.linspace(0.05,0.2,100)
d = np.array([0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 10])
A_H = 1
N = 100
hc = 0.5
#plot_fold(5,1,D,d,hc,A_H)

# Probamos
#a = 1
#A_H = 1
#d = 0.1
#r = 0.03
#hc = 0.5
#fold = 5
#N = 1
#n = 1
#ico = Polyhedron.icosahedron(fold = 5)
#l = 100
#D = np.linspace(0.05, 0.2, l)
#
#F_1 = np.empty(l)
#F   = np.empty(l)
#F_punt = np.empty(l)
#F_tot = np.empty(l)
#
#for i in range(l):
#
#    F[i]      = vdW_face(fold,a,hc,D[i],d,A_H)[2]
#    F_punt[i] = point_vdw(ico, r, D[i], A_H, N, n , hc, True)[2]
#    F_1[i]    = vdW_face_r_good(fold,a,r,hc,D[i],d,A_H)[2]
#
#F_tot = F + F_punt
##print(F_punt)
#
#plt.plot(D,F,'.')
#plt.plot(D,F_1,'.')
#plt.xlabel('Minimum distance (D)')
#plt.ylabel('Fvdw')
#plt.title('Fvdw with one protein at face center (2 methods)')
#plt.grid()
#plt.legend(['Complete', 'Old'])
#plt.show()

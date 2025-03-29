from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import Derj_IcoIco as derj



PHI = (1 + np.sqrt(5)) / 2


# Function that creates the list of lists used for the rotations
def angles(Fold_in, Fold_out, N):
    
    phi = (1+np.sqrt(5))/2
    out = []

    if Fold_in == 5:
        
        if Fold_out == 3:
            ang_x = np.linspace(0,38.4*np.pi/180,N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return
        
        elif Fold_out == 2:
            ang_x = np.linspace(0,-np.arctan(phi/(phi+1)),N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return out

        
    if Fold_in == 3:
        
        if Fold_out == 5:
            ang_x = np.linspace(0,-38.4*np.pi/180,N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return out
        
        elif Fold_out == 2:
            ang_x = np.linspace(0,-np.arctan(phi/(phi-1)),N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return out

    if Fold_in == 2:
        
        if Fold_out == 5:
            ang_x = np.linspace(0,np.arctan(phi/(phi+1)),N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return out
        
        elif Fold_out == 3:
            ang_x = np.linspace(0,np.arctan(phi/(phi-1)),N)
            for i in range(N):
                out.append([ang_x[i],0,0])
            return out


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


# Function that, for a fixed distance D, plane depth d and considered height hc, calculates and plots the Fvdw along the 
# rotation path from one fold to another
def F_rot_comp(Fold1, Fold2, Fold3, a, D, d, hc, A_H, N):

    rotation1 = angles(Fold1, Fold2, N)     # rot from fold1 to fold2
    rotation2 = angles(Fold2, Fold3, N)     # rot from fold2 to fold3

    rotation = rotation1 + rotation2        # just for plotting




def F_rot(Fold1, Fold2, a, D, d, hc, A_H, N):


    F = []
    rotation = angles(Fold1, Fold2, N)

    for i in range(N):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=Fold1, rot_extra=rotation[i])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        areas, normals = Polyhedron.area_normal(ico, hc+correction, a)

        areas = np.array(areas)
        normals = np.array(normals)
        F_temp = np.empty([len(areas),3])

        for j in range(len(areas)):
            F_module = abs(derjaguin_polygon(areas[j], D, d, A_H))
            proj = np.array(projection(normals[j])[2,:])
            F_temp[j,:] = np.array(F_module*proj)

        F_add = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
        F.append(F_add)
    
    return np.array(F)

def F_rot1(Fold1, a, D, d, hc, A_H, N):


    F = []
    rotation = angles1_D(N)

    for i in range(N):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=Fold1, rot_extra=rotation[i])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        areas, normals = Polyhedron.area_normal(ico, hc+correction, a)

        areas = np.array(areas)
        normals = np.array(normals)
        F_temp = np.empty([len(areas),3])

        for j in range(len(areas)):
            F_module = abs(derjaguin_polygon(areas[j], D, d, A_H))
            proj = np.array(projection(normals[j])[2,:])
            F_temp[j,:] = np.array(F_module*proj)

        F_add = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
        F.append(F_add)
    
    return np.array(F)



# Function to plot the rotation from 2-Fold to 5-fold
def rot2_3(a, D, d, hc, A_H, N):

    F = F_rot(2,3,a,D,d,hc,A_H,N)
    rot_max = angles(2,3,N)[-1][0]
    rot = np.linspace(0, rot_max, N)

    F_x = F[:,0]
    F_y = F[:,1]
    F_z = F[:,2]

    #plt.plot(rot, F_x,'.')
    #plt.plot(rot, F_y,'.')
    plt.plot(rot, F_z,'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-3 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['Fz'])
    plt.grid()
    plt.show()

def rot3_5(a, D, d, hc, A_H, N):

    F = F_rot(3,5,a,D,d,hc,A_H,N)
    rot_max = angles(3,5,N)[-1][0]
    rot = np.linspace(0, rot_max, N)

    F_x = F[:,0]
    F_y = F[:,1]
    F_z = F[:,2]

    #plt.plot(rot, F_x,'.')
    #plt.plot(rot, F_y,'.')
    plt.plot(rot, F_z,'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-3 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['Fz'])
    plt.grid()
    plt.show()


# Function to plot the rotation from 2-Fold to 5-fold
def rot2_5(a, D, d, hc, A_H, N):

    F = F_rot(2,5,a,D,d,hc,A_H,N)
    rot_max = angles(2,5,N)[-1][0]
    rot = np.linspace(0, rot_max, N)

    F_x = F[:,0]
    F_y = F[:,1]
    F_z = F[:,2]

    #plt.plot(rot, F_x,'.')
    #plt.plot(rot, F_y,'.')
    plt.plot(rot, F_z,'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-5 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['Fz'])
    plt.grid()
    plt.show()


# Rotation from 2-fold to 3-fold, to 5-fold and back to 2-fold
def rot_full(a, D, d, hc, A_H, N):

    F1 = F_rot(2,3,a,D,d,hc,A_H,N)[:,2]
    F2 = F_rot(3,5,a,D,d,hc,A_H,N)[:,2]
    F3 = F_rot(5,2,a,D,d,hc,A_H,N)[:,2]

    F_ = np.append(F1,F2)
    F = np.append(F_,F3)
    
    rot = np.linspace(0,3,3*N)

    plt.plot(rot, F,'.')
    plt.axvline(0,color='r')
    plt.axvline(1,color='r')
    plt.axvline(2,color='r')
    plt.axvline(3,color='r')
    plt.title('Forces in axis (x,y,z) vs rotation 2-5 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['F','2-fold','3-fold','5-fold','2-fold'])
    plt.grid()
    plt.savefig('rotaciones.png')




def rot2_5_hc(a, D, d, hc, A_H, N):
    
    F = np.empty([5,N])
    
    for i in range(5):
        F[i,:] = F_rot(2,5,a,D,d,hc[i],A_H,N)[:,2]

    rot_max = angles(2,5,N)[-1][0]
    rot = np.linspace(0, rot_max, N)

    for i in range(5):
        plt.plot(rot, F[i,:],'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-5 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['hc = 0.1','hc = 0.25','hc = 0.5','hc = 0.75','hc = 1'])
    plt.grid()
    plt.show()

def rot2_3_hc(a, D, d, hc, A_H, N):
    
    F = np.empty([5,N])
    
    for i in range(5):
        F[i,:] = F_rot(2,3,a,D,d,hc[i],A_H,N)[:,2]

    rot_max = angles(2,3,N)[-1][0]
    rot = np.linspace(0, rot_max, N)

    for i in range(5):
        plt.plot(rot, F[i,:],'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-3 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['hc = 0.1','hc = 0.25','hc = 0.5','hc = 0.75','hc = 1'])
    plt.grid()
    plt.show()


def angles1_D(N):
    out = []
    ang_x = np.linspace(0,np.pi/2,N)
    for i in range(N):
        out.append([ang_x[i],0,0])

    return out



def rot_1D(a, D, d, hc, A_H, N):

    F = np.empty([5,N])
    
    for i in range(5):
        F[i,:] = F_rot1(5,a,D,d,hc[i],A_H,N)[:,2]

    rot = angles1_D(N)
    
    for i in range(5):
        plt.plot(rot, F[i,:],'.')
    plt.title('Forces in axis (x,y,z) vs rotation 2-5 fold')
    plt.xlabel('Rotation angle (in axis x)')
    plt.ylabel('Van der Waals forces')
    plt.legend(['hc = 0.1','hc = 0.25','hc = 0.5','hc = 0.75','hc = 1'])
    plt.grid()
    plt.show()




def U_5_3_2(a, D, d, hc, A_H, N):

    rot5_3 = np.linspace(0, 38.4*np.pi/180, N)
    rot3_2 = np.linspace(38.4*np.pi/180, (38.4+18.33)*np.pi/180, N) 
    rotation = np.append(rot5_3,rot3_2)

    F = []

    for i in range(len(rotation)):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=5, rot_extra=[rotation[i],0,0])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        areas, normals = Polyhedron.area_normal(ico, hc+correction, a)

        areas = np.array(areas)
        normals = np.array(normals)
        F_temp = np.empty([len(areas),3])

        for j in range(len(areas)):
            F_module = abs(derjaguin_polygon(areas[j], D, d, A_H))
            proj = np.array(projection(normals[j])[2,:])
            F_temp[j,:] = np.array(F_module*proj)

        F_add = np.array([np.sum(F_temp[:,0]), np.sum(F_temp[:,1]), np.sum(F_temp[:,2])])
        F.append(F_add)
    F = np.array(F)
    #U = simpson(F[:,2], rotation)
    #print(U)

    return F[:,2]

def rot5_3_2(a,D,d,hc,A_H,N,one):

    if one == True:
        U = U_5_3_2(a,D,d,hc,A_H,N)
        rot5_3 = np.linspace(0, 38.4, N)
        rot3_2 = np.linspace(38.4, 38.4+18.33, N) 
        rotation = np.append(rot5_3,rot3_2)
        plt.plot(rotation,U,'.')
        plt.title('Fuerza de vdw en la rotación 5-3-2 (eje X)')
        plt.xlabel('Ángulo de rotación (º)')
        plt.ylabel('Fuerza de van der Waals ($A_H$/a)')
        plt.axvline(x=0, color = 'b')
        plt.axvline(x=37.4, color = 'tab:orange')
        plt.axvline(x=38.4+18.33, color = 'g')
        plt.legend(['Fuerza','5-fold','3-fold','2-fold'])
        plt.grid()
        #plt.savefig('Rot_1.png')
        plt.show()

    elif one ==2:
        U = U_5_3_2(a,D,d,hc,A_H,N)
        
        rot5_3 = np.linspace(0, 38.4, N)
        rot3_2 = np.linspace(38.4, 38.4+18.33, N) 
        rotation = np.append(rot5_3,rot3_2)

        energ = simpson(U,rotation)
    

        plt.plot(rotation,U,'.')
        plt.title('Fuerza de vdw en la rotación 5-3-2 (eje X)')
        plt.xlabel('Ángulo de rotación (º)')
        plt.ylabel('Fuerza de van der Waals ($A_H$/a)')
        plt.axvline(x=0, color = 'b')
        plt.axvline(x=37.4, color = 'tab:orange')
        plt.axvline(x=38.4+18.33, color = 'g')
        plt.legend(['Fuerza','5-fold','3-fold','2-fold'])
        plt.grid()
        #plt.savefig('Rot_1.png')
        plt.show()

    elif one == 3:
        rot5_3 = np.linspace(0, 38.4, N)
        rot3_2 = np.linspace(38.4, 38.4+18.33, N) 
        rotation = np.append(rot5_3,rot3_2)

        U = np.empty([len(hc),len(rotation)])
        for i in range(len(hc)):
            F = U_5_3_2(a,D,d,hc[i],A_H,N)
            energ = np.array(simpson(F,rotation))
            U[i,:] = energ


        #for i in range(len(hc)):
        #    plt.plot(rotation, U[i,:],'.')
        plt.plot(rotation,U[1,:],'.')
        plt.title('Energías de vdw en la rotación 5-3-2 (eje X)')
        plt.xlabel('Ángulo de rotación (º)')
        plt.ylabel('Energía de van der Waals ($A_H$/a)')
        plt.axvline(x=0, color = 'b')
        plt.axvline(x=37.4, color = 'tab:orange')
        plt.axvline(x=38.4+18.33, color = 'g')
        plt.legend(['Energía','5-fold','3-fold','2-fold'])
        plt.legend(['hc = 0.7','hc = 0.768','hc = 0.85'])
        plt.grid()
        plt.show()
        #plt.savefig('Rot_cool.png')
        
    else:
        rot5_3 = np.linspace(0, 38.4, N)
        rot3_2 = np.linspace(38.4, 38.4+18.33, N) 
        rotation = np.append(rot5_3,rot3_2)

        U = np.empty([len(hc),len(rotation)])
        for i in range(len(hc)):
            U[i,:] = U_5_3_2(a,D,d,hc[i],A_H,N)

        
        for i in range(len(hc)):
            plt.plot(rotation, U[i,:],'.')
        plt.title('Fuerzas de vdw en la rotación 5-3-2 (eje X)')
        plt.xlabel('Ángulo de rotación (º)')
        plt.ylabel('Fuerza de van der Waals ($A_H$/a)')
        plt.axvline(x=0, color = 'b')
        plt.axvline(x=37.4, color = 'tab:orange')
        plt.axvline(x=38.4+18.33, color = 'g')
        plt.legend(['Fuerza','5-fold','3-fold','2-fold'])
        plt.legend(['hc = 0.7','hc = 0.768','hc = 0.85'])
        plt.grid()
        plt.show()
        #plt.savefig('Rot_cool.png')

    return
a = 1
D = 0.06
d = 0.1
#hc = [0.1, 0.25, 0.5, 0.75, 1]
hc = 0.25
#hc = [0.7,0.768,0.85]
A_H = 1
N = 100

#rot5_3_2(a,D,d,hc,A_H,N,True)



#ico = Polyhedron.icosahedron(fold = 2, rot_extra=[np.arctan(PHI/(PHI+1)),0,0])
#ico.plot(show_labels=True)
#rot3_5(a,D,d,0.5,A_H,N)
#rot_1D(a,D,d,hc,A_H,N)
#rot_full(a,D,d,0.5,A_H,N)


# Function that calculates the average distance between the particles forming the virus and a plane situated at a distance D
def dist_med(a, D):
    
    # F:  fold
    # D:  minimum distance between icosahedron and plane
    # a:  length of icosahedron edge
    # k:  wave number of surface
    # hm: max height of surface
    # hc: maximum height considered (if it's more that the one allowed by the script it becomes that, the allowed heights are
    # described in hc_max)


    R = a*np.sin(2*np.pi/5)
    M = 100
    hc = 0.4*a

    rot5_3 = np.linspace(0, 38.4*np.pi/180, N)
    rot3_2 = np.linspace(38.4*np.pi/180, (38.4+18.33)*np.pi/180, N) 
    rotation = np.append(rot5_3,rot3_2)

    d_m = []
    for i in range(len(rotation)):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=5, rot_extra=[rotation[i],0,0])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        ico2 = Polyhedron.icosahedron(edge_length=a, fold = 5, rot_extra=[rotation[i],0,0], trans_extra=[0,0,-correction])

        # Defino los puntos
        puntos = ico2.get_points()

        # Defino las caras
        faces = fc.triangulos(puntos)

        # Creo los puntos sobre las caras
        puntos_face = []
        for i in range(len(faces)):
            puntos_face.extend(fc.planos(faces[i],M))

        # Selecciono los puntos que me interesan (los que están por debajo de la altura de corte)
        puntos_temp = fc.CorteZ(puntos_face,hc)

        puntos = fc.Transf(puntos_temp,[0,0,D],[0,0,0])
        supV = fc.SupV(puntos,R,M)

        d_m.append(np.average(supV))

    return np.array(d_m)

# D no conservada
def dist_med_D(a, D):
    
    # F:  fold
    # D:  minimum distance between icosahedron and plane
    # a:  length of icosahedron edge
    # k:  wave number of surface
    # hm: max height of surface
    # hc: maximum height considered (if it's more that the one allowed by the script it becomes that, the allowed heights are
    # described in hc_max)


    R = a*np.sin(2*np.pi/5)
    M = 100
    hc = 0.4*a

    rot5_3 = np.linspace(0, 38.4*np.pi/180, N)
    rot3_2 = np.linspace(38.4*np.pi/180, (38.4+18.33)*np.pi/180, N) 
    rotation = np.append(rot5_3,rot3_2)

    d_m = []
    for i in range(len(rotation)):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=5, rot_extra=[rotation[i],0,0], trans_extra=[0,0,D])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        ico2 = Polyhedron.icosahedron(edge_length=a, fold = 5, rot_extra=[rotation[i],0,0], trans_extra=[0,0,-correction])

        # Defino los puntos
        puntos = ico2.get_points()

        # Defino las caras
        faces = fc.triangulos(puntos)

        # Creo los puntos sobre las caras
        puntos_face = []
        for i in range(len(faces)):
            puntos_face.extend(fc.planos(faces[i],M))

        # Selecciono los puntos que me interesan (los que están por debajo de la altura de corte)
        puntos_temp = fc.CorteZ(puntos_face,hc)

        puntos = fc.Transf(puntos_temp,[0,0,D],[0,0,0])
        supV = fc.SupV(puntos,R,M)

        d_m.append(np.average(supV))

    return np.array(d_m)
    
def plot_dist(d_m):

    N = 100
    rot5_3 = np.linspace(0, 38.4, N)
    rot3_2 = np.linspace(38.4, (38.4+18.33), N) 
    rotation = np.append(rot5_3,rot3_2)

    plt.plot(rotation,d_m,'.')
    plt.title('Average distance between virus and planar surface')
    plt.xlabel('Rotation angle (º)')
    plt.ylabel('Average distance (a)')
    plt.axvline(x=0, color = 'b')
    plt.axvline(x=37.4, color = 'tab:orange')
    plt.axvline(x=38.4+18.33, color = 'g')
    plt.legend(['Average distance','5-fold','3-fold','2-fold'])
    plt.grid()
    plt.show()

    return


# Area de interaccion
def Area(a, hc, N):

    rot5_3 = np.linspace(0, 38.4*np.pi/180, N)
    rot3_2 = np.linspace(38.4*np.pi/180, (38.4+18.33)*np.pi/180, N) 
    rotation = np.append(rot5_3,rot3_2)

    area = []

    for i in range(len(rotation)):

        # create icosahedron
        ico   = Polyhedron.icosahedron(edge_length=a ,fold=5, rot_extra=[rotation[i],0,0])  
        
        # calculate correction so that hc starts from min_point of icosahedron
        correction = Polyhedron.z_min(ico)   
                  
        # now we just calculate the force
        areas, normals = Polyhedron.area_normal(ico, hc+correction, a)

        areas = np.array(areas)
        area.append(sum(areas))

    return area

def plot_area(a,hc,N):
    
    area = Area(a, hc, N)
    rot5_3 = np.linspace(0, 38.4, N)
    rot3_2 = np.linspace(38.4, (38.4+18.33), N) 
    rotation = np.append(rot5_3,rot3_2)
   

    plt.plot(rotation,area,'.')
    plt.title('Interacting area of the virus vs rotation angles')
    plt.xlabel('Rotation angle (º)')
    plt.ylabel('Interacting area (a²)')
    plt.axvline(x=0, color = 'b')
    plt.axvline(x=37.4, color = 'tab:orange')
    plt.axvline(x=38.4+18.33, color = 'g')
    plt.legend(['Area','5-fold','3-fold','2-fold'])
    plt.grid()
    plt.show()


    


a = 1
D = 0.1
hc = 0.4
N = 100
#plot_area(a,hc,100)






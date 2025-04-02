from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt

# The objective of this set of functions is plotting the van der Waals force between an icosahedron in different symmetries and
# a plane with respect to the minimum distance between them

# Function that calculates the area of the cross section of the icosahedron at certain height 
def area_trans(ico, z):

    # ico: icosahedron created with the desired characteristics 
    # z: height at which the cross section is considered

    surface = ico.intersect(z)
    #print(surface)
    return surface.area


def a_cross_sphere(R, z):

    # R: radious of sphere 
    # z: height at which the cross section is considered
    area = np.pi*(R*R-(R-z)**2)

    return area

# Function that calculates the derivative of a given function using the Finite Difference Method (with error O(z^4))
def deriv(A_z, z):

    # A_z: list containing the values for the cross sections 
    # z: numpy array containing the corresponding z of the cross sections
    
    A_z = np.array(A_z)
    dA_z = np.empty(len(A_z))
    h = z[1] - z[0]

    # previous derivative on the boundary
    # dA_z[0]  = (A_z[1]  - A_z[0])/h
    # dA_z[1]  = (A_z[2]  - A_z[0])/(2*h)
    # dA_z[-2] = (A_z[-1] - A_z[-3])/(2*h)
    # dA_z[-1] = (A_z[-1] - A_z[-2])/h

    # First and last two terms with fourth order schemes
    dA_z[0]  = ((-1/4)*A_z[4] + (4/3)*A_z[3] + (-3)*A_z[2] + (4)*A_z[1] + (-25/12)*A_z[0])/h 
    dA_z[1]  = ((-1/4)*A_z[5] + (4/3)*A_z[4] + (-3)*A_z[3] + (4)*A_z[2] + (-25/12)*A_z[1])/h 
    dA_z[-2]  = ((1/4)*A_z[-6] + (4/3)*A_z[-5] + (3)*A_z[-4] + (4)*A_z[-3] + (25/12)*A_z[-2])/h 
    dA_z[-1]  = ((1/4)*A_z[-5] + (4/3)*A_z[-4] + (3)*A_z[-3] + (4)*A_z[-2] + (25/12)*A_z[-1])/h 

    # Rest of the terms
    for i in np.arange(2,len(A_z)-2):
        dA_z[i] = (A_z[i-2] - 8*A_z[i-1] + 8*A_z[i+1] - A_z[i+2])/(12*h)

    return dA_z


# Function that numerically calculates the Darjeguin approximation to the van der Waals force between an icosahedron 
# (in 2, 3 or 5-fold symmetry) with a minimum distance of D = min_height, using Scipy's integrated Simpson's rule for 
# the integration
def darjeguin(Fold, h_cutoff, edge_length, min_height, A_H, N, double):

    # Fold: symmetry being used, must be an integer (2, 3 or 5).
    # h_cutoff: maximum height being considered for the interaction, float (must be between min_height and the maximum height
    # of the icosahedron). Note: h_cutoff is defined with 0 being the surface of the plane, contrary to the convention used 
    # in the other scripts.
    # edge_length: length of the icosahedron's edges, image.pngfloat.
    # min_height: minimum distance between the icosahedron and the plane (also called D), float.
    # A_H: Hamaker constant of the system, float
    # N: number of intervals (cross sections) considered between min_height and h_cutoff, integer 
    # double: indica si se considera la figura duplicada, o un plano (True - figura, False - plano)

    # First step: creating the icosahedron and positioning it at the correct height
    
    # Extra rotation in order to allow for a smooth derivative for the 3-fold
    if Fold == 3:
        rot = 0.017
        ico_temp = Polyhedron.icosahedron(edge_length = edge_length, fold = Fold, trans_extra=[0,0,0], rot_extra=[rot,0,0])
        z_extra = Polyhedron.z_min(ico_temp)
        ico = Polyhedron.icosahedron(edge_length = edge_length, fold = Fold, trans_extra=[0,0,min_height-z_extra], rot_extra=[rot,0,0])

    else:
        ico = Polyhedron.icosahedron(edge_length = edge_length, fold = Fold, trans_extra=[0,0,min_height], rot_extra=[0,0,0])
    

    # Second step: calculating A(z) for every z included in the interval (min_height, hc)
    dz = (h_cutoff-min_height)/(N-1)
    z = np.linspace(min_height, h_cutoff, N)
    z_ant = np.array([min_height-dz])

    # Mimics the fact that for a doubled figure, the distance x required to reach a certain area is double that which is 
    # needed in a figure vs plane setting
    if double == True:
        z = z/2

    A_z = [0]
    for i in z:
        A_z.append(area_trans(ico, i))
    
    if double == True:
        z = 2*z

    z = np.append(z_ant,z)

    # Third sep: deriving A(z) and adding the other term in the Darjeguin approximation
    dA_z = deriv(A_z, z)
    f_z = -(A_H/(6*np.pi*z*z*z))*dA_z
    
    #C1 = np.sqrt((5-np.sqrt(5))/10)
    #A_teo = 0.25*np.sqrt(5*(5+2*np.sqrt(5)))*((z-min_height)/C1)**2
    #plt.plot(z,A_z)
    #plt.plot(z,A_teo,'.')
    #plt.plot(z,dA_z)
    #plt.legend(['A_z','A_teo','dA_z'])
    #plt.show()

    # Fourth step: computing the integral of f_z from min_height to h_cutoff
    F_D = simpson(f_z, z)
    
    return F_D

def derjaguin_sphere(R, min_height, h_cutoff, A_H, N, double):

    # h_cutoff: maximum height being considered for the interaction, float (must be between min_height and the maximum height
    # of the sphere). Note: h_cutoff is defined with 0 being the surface of the plane, contrary to the convention used 
    # in the other scripts.
    # min_height: minimum distance between the sphere and the plane (also called D), float.
    # A_H: Hamaker constant of the system, float
    # N: number of intervals (cross sections) considered between min_height and h_cutoff, integer 

    # First step: calculating A(z) for every z included in the interval (min_height, max_height) for both approximations
    z = np.linspace(min_height, h_cutoff, N) - min_height

    if double == True:
        z = z/2

    A_z_cross = a_cross_sphere(R,z)

    if double == True:
        z = 2*z

    # Second sep: deriving A(z) and adding the other term in the Darjeguin approximation
    dA_z_cross = deriv(A_z_cross, z)

    # Third step: undoing the change to z and adding the term VA(z) to the function
    x = z + min_height
    f_x_cross = - (A_H/(6*np.pi*x*x*x))*dA_z_cross

    # Fourth step: computing the integral of f_z_() from min_height to h_cutoff
    F_D_cross = simpson(f_x_cross, x)
    
    return F_D_cross


def main(F,D):
    N = 100 
    a = 1
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)
    coso2 = np.linspace(0.01,0.2,N)

    # 5-fold with analytical solution
    if F == 6:
        h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))
        k = coso2 + h
        F_D = np.empty(len(coso2))
        for i in range(len(coso2)):
            F_D[i] = darjeguin(Fold=5, h_cutoff=coso2[i] + h, edge_length=a, min_height=coso2[i], A_H=1, N=N, double=False)
        C1 = np.sqrt((5-np.sqrt(5))/10)
        C2 = 0.5*np.sqrt(5*(5+2*np.sqrt(5)))
        F_teo = - C2/(12*np.pi*C1**2)*(1/coso2-(2*k-coso2)/k**2)

        plt.plot(coso2,F_D,'.')
        plt.plot(coso2,F_teo,'.')
        plt.legend(['Calculado','Analítico'])
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de vand der Waals (U.A.)')
        plt.title('Fvdw entre pirámide pentagonal y plano')
        plt.grid()
        plt.show()

    # Relative error between analytical and numerical results in the previous case
    elif F == 0:
        h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))
        k = coso2 + h
        F_D = np.empty(len(coso2))
        for i in range(len(coso2)):
            F_D[i] = darjeguin(Fold=5, h_cutoff=coso2[i] + h, edge_length=a, min_height=coso2[i], A_H=1, N=N, double=False)
        C1 = np.sqrt((5-np.sqrt(5))/10)
        C2 = 0.5*np.sqrt(5*(5+2*np.sqrt(5)))
        F_teo = - C2/(12*np.pi*C1**2)*(1/coso2-(2*k-coso2)/k**2)
        
        F_err = np.abs(F_D - F_teo)/np.abs(F_D)
        plt.plot(coso2,F_err,'.')
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Relative error')
        plt.title('Relative error between analytical and numerical results')
        plt.show()

    # All folds + sphere
    elif F == 1:
        F_D2 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D2[i] = darjeguin(Fold=2, h_cutoff=coso[i]+hc[2], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
        F_D3 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D3[i] = darjeguin(Fold=3, h_cutoff=coso[i]+hc[1], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
        F_D5 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D5[i] = darjeguin(Fold=5, h_cutoff=coso[i]+hc[0], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
        F_DS = np.empty(len(coso))
        for i in range(len(coso)):
            F_DS[i] = derjaguin_sphere(R=R, min_height=coso[i], h_cutoff = coso[i] + 2*R, A_H=1, N=N, double=False)
        
        coso = np.log(coso)
        F_D2 = -np.log(-F_D2)
        F_D3 = -np.log(-F_D3)
        F_D5 = -np.log(-F_D5)

        plt.plot(coso,F_D2,'.')
        plt.plot(coso,F_D3,'.')
        plt.plot(coso,F_D5,'.')
        plt.plot(coso,np.log(F_DS))
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de interacción vdW')
        plt.title('Fuerza de van der Waals entre cada simetría y el plano')
        plt.legend(['2-fold','3-fold','5-fold','sphere'])
        plt.show()

    # All folds
    elif F == 7:
        F_D2 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D2[i] = darjeguin(Fold=2, h_cutoff=coso[i]+hc[2], edge_length=a, min_height=coso[i], A_H=1, N=N)
        F_D3 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D3[i] = darjeguin(Fold=3, h_cutoff=coso[i]+hc[1], edge_length=a, min_height=coso[i], A_H=1, N=N)
        F_D5 = np.empty(len(coso))
        for i in range(len(coso)):
            F_D5[i] = darjeguin(Fold=5, h_cutoff=coso[i]+hc[0], edge_length=a, min_height=coso[i], A_H=1, N=N)

        plt.plot(coso,F_D2,'.')
        plt.plot(coso,F_D3,'.')
        plt.plot(coso,F_D5,'.')
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de interacción vdW')
        plt.title('Fuerza de van der Waals entre cada simetría y el plano')
        plt.legend(['2-fold','3-fold','5-fold'])
        plt.show()

    # Sphere-Plane vs Sphere-Sphere
    elif F == 8:
        F_S = np.empty(len(coso))
        for i in range(len(coso)):
            F_S[i] = derjaguin_sphere(R=R, min_height=coso[i], h_cutoff = coso[i] + 2*R, A_H=1, N=N, double=False)
        F_2S = np.empty(len(coso))
        for i in range(len(coso)):
            F_2S[i] = derjaguin_sphere(R=R, min_height=coso[i], h_cutoff = coso[i] + 2*R, A_H=1, N=N, double=True)
        
        plt.plot(coso,F_S,'.')
        plt.plot(coso,F_2S,'.')
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de interacción vdW')
        plt.title('Fuerza de van der Waals Esfera-Plano y Esfera-Esfera')
        plt.legend(['Esfera-Plano','Esfera-Esfera'])
        plt.show()

    # Fvdw vs hc for each simmetry + sphere (for D = 0.025)
    elif F == 10:
        #D = 0.025
        N = 300
        hc2 = np.linspace(0.01,hc[2],100)
        hc3 = np.linspace(0.01,hc[1],100)
        hc5 = np.linspace(0.01,hc[0],100)
        #hcS = np.linspace(0.01,2*R,100)

        F_2 = np.empty(len(hc2))
        F_3 = np.empty(len(hc2))
        F_5 = np.empty(len(hc2))
        #F_S = np.empty(len(hc2))

        for i in range(len(hc2)):
            F_2[i] = darjeguin(Fold=2, h_cutoff=D+hc2[i], edge_length=a, min_height=D, A_H=1, N=N, double=False)
            F_3[i] = darjeguin(Fold=3, h_cutoff=D+hc3[i], edge_length=a, min_height=D, A_H=1, N=N, double=False)
            F_5[i] = darjeguin(Fold=5, h_cutoff=D+hc5[i], edge_length=a, min_height=D, A_H=1, N=N, double=False)
            #F_S[i] = derjaguin_sphere(R=R, min_height=D, h_cutoff =D+hcS, A_H=1, N=N, double=True)
        
        plt.plot(hc5,F_5,'.')
        plt.plot(hc3,F_3,'.')
        plt.plot(hc2,F_2,'.')
        #plt.plot(hcS,F_S,'.')
        plt.xlabel('Maximum height considered in figure')
        plt.ylabel('Van der Waals force plane-figure')
        plt.title('VdW force vs considered height')
        plt.legend(['Ico 5-Fold','Ico 3-Fold','Ico 2-Fold'])
        plt.show()

    elif F == 11:
        Dp = np.linspace(0.1,1,10)
        N = 100

        hc5 = np.linspace(0.01,hc[0],100)
        F_5 = np.empty([len(Dp),len(hc5)])

        for j in range(len(Dp)):
            for i in range(len(hc5)):
                F_5[j,i] = darjeguin(Fold=5, h_cutoff=Dp[j]+hc5[i], edge_length=a, min_height=Dp[j], A_H=1, N=N, double=False)

        for i in range(len(Dp)):
            plt.plot(hc5,F_5[i,:],'.')
        plt.xlabel('Maximum height considered in figure')
        plt.ylabel('Van der Waals force plane-figure')
        plt.title('VdW force vs considered height for 5-fold')
        plt.legend(['D = 0.1','D = 0.2','D = 0.3','D = 0.4','D = 0.5','D = 0.6','D = 0.7','D = 0.8','D = 0.9','D = 1.0'],loc='best')
        plt.grid()
        plt.show()
    
    elif F == 12:
        Dp = np.linspace(0.1,1,10)
        N = 100

        hc3 = np.linspace(0.01,hc[1],100)
        F_3 = np.empty([len(Dp),len(hc3)])

        for j in range(len(Dp)):
            for i in range(len(hc3)):
                F_3[j,i] = darjeguin(Fold=3, h_cutoff=Dp[j]+hc3[i], edge_length=a, min_height=Dp[j], A_H=1, N=N, double=False)

        for i in range(len(Dp)):
            plt.plot(hc3,F_3[i,:],'.')
        plt.xlabel('Maximum height considered in figure')
        plt.ylabel('Van der Waals force plane-figure')
        plt.title('VdW force vs considered height for 3-fold')
        plt.legend(['D = 0.1','D = 0.2','D = 0.3','D = 0.4','D = 0.5','D = 0.6','D = 0.7','D = 0.8','D = 0.9','D = 1.0'],loc='best')
        plt.grid()
        plt.show()
    
    elif F == 13:
        Dp = np.linspace(0.1,1,10)
        N = 100

        hc2 = np.linspace(0.01,hc[2],100)
        F_2 = np.empty([len(Dp),len(hc2)])

        for j in range(len(Dp)):
            for i in range(len(hc2)):
                F_2[j,i] = darjeguin(Fold=2, h_cutoff=Dp[j]+hc2[i], edge_length=a, min_height=Dp[j], A_H=1, N=N, double=False)

        for i in range(len(Dp)):
            plt.plot(hc2,F_2[i,:],'.')
        plt.xlabel('Maximum height considered in figure')
        plt.ylabel('Van der Waals force plane-figure')
        plt.title('VdW force vs considered height, for 2-fold')
        plt.legend(['D = 0.1','D = 0.2','D = 0.3','D = 0.4','D = 0.5','D = 0.6','D = 0.7','D = 0.8','D = 0.9','D = 1.0'],loc='best')
        plt.grid()
        plt.show()

    else:
        F_D = np.empty(len(coso))
        hc_in = 0*(F==5) + 1*(F==3) + 2*(F==2)
        for i in range(len(coso)):
            F_D[i] = darjeguin(Fold=F, h_cutoff=coso[i]+hc[hc_in], edge_length=a, min_height=coso[i], A_H=1, N=N)
        plt.plot(coso,F_D)
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de interacción vdW')
        plt.title('Fuerza de van der Waals entre simetría ' + str(F) +'-fold y el plano')
        plt.show()
    

    
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

# Function that calculates the derivative of a given function using the Finite Difference Method (with error O(z^4))
def deriv(A_z, z):

    # A_z: list containing the values for the cross sections 
    # z: numpy array containing the corresponding z of the cross sections

    
    A_z = np.array(A_z)
    dA_z = np.empty(len(A_z))
    h = z[1] - z[0]

    # First and last two terms
    dA_z[0]  = (A_z[1]  - A_z[0])/h
    dA_z[1]  = (A_z[2]  - A_z[0])/(2*h)
    dA_z[-2] = (A_z[-1] - A_z[-3])/(2*h)
    dA_z[-1] = (A_z[-1] - A_z[-2])/h

    # Rest of the terms
    for i in np.arange(2,len(A_z)-2):
        dA_z[i] = (A_z[i-2] - 8*A_z[i-1] + 8*A_z[i+1] - A_z[i+2])/(12*h)

    return dA_z


# Function that numerically calculates the Darjeguin approximation to the van der Waals force between an icosahedron 
# (in 2, 3 or 5-fold symmetry) with a minimum distance of D = min_height, using Scipy's integrated Simpson's rule for 
# the integration
def darjeguin(Fold, h_cutoff, edge_length, min_height, A_H, N):

    # Fold: symmetry being used, must be an integer (2, 3 or 5).
    # h_cutoff: maximum height being considered for the interaction, float (must be between min_height and the maximum height
    # of the icosahedron). Note: h_cutoff is defined with 0 being the surface of the plane, contrary to the convention used 
    # in the other scripts.
    # edge_length: length of the icosahedron's edges, float.
    # min_height: minimum distance between the icosahedron and the plane (also called D), float.
    # A_H: Hamaker constant of the system, float
    # N: number of intervals (cross sections) considered between min_height and h_cutoff, integer 


    # First step: creating the icosahedron and positioning it at the correct height
    ico = Polyhedron.icosahedron(edge_length = edge_length, fold = Fold, trans_extra=[0,0,min_height])
    #ico.plot()


    # Second step: calculating A(z) for every z included in the interval (min_height, hc)
    z = np.linspace(min_height, h_cutoff, N)
    A_z = []
    for i in z:
        A_z.append(area_trans(ico, i))
    
    
    
    
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

def main(F):
    N = 100 
    a = 1
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.01, 0.2, N)

    # 5-fold with analytical solution
    if F == 6:
        h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))
        k = coso + h
        F_D = np.empty(len(coso))
        for i in range(len(coso)):
            F_D[i] = darjeguin(Fold=5, h_cutoff=coso[i] + h, edge_length=a, min_height=coso[i], A_H=1, N=N)
        C1 = np.sqrt((5-np.sqrt(5))/10)
        C2 = 0.5*np.sqrt(5*(5+2*np.sqrt(5)))
        F_teo = - C2/(12*np.pi*C1**2)*(1/coso-(2*k-coso)/k**2)
        plt.plot(coso,F_D,'.')
        plt.plot(coso,F_teo,'.')
        plt.legend(['calculado','analitico'])
        plt.xlabel('Distancia mínima D (a)')
        plt.ylabel('Fuerza de interacción vdW')
        plt.title('Fuerza de van der Waals entre simetría 5 fold (pirámide) y el plano')
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






















































































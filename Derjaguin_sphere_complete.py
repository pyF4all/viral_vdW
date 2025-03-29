from shapely.geometry import Point, Polygon
import numpy as np
from scipy.integrate import simpson
from scipy.spatial.transform import Rotation as R
import funciones_final as fc
from polyhedron_final import  Polyhedron
import math as math
import matplotlib.pyplot as plt

# The objective of this set of functions is plotting the van der Waals force between a semisphere and
# a plane with respect to the minimum distance between them, and comparing the results obtained between the use of 
# Darjeguin's aproximation and a modification of it where we consider surface area instead of cross section area

# Function that calculates the area of the cross section of the sphere at certain height 
def a_cross_sphere(R, z):

    # R: radious of sphere 
    # z: height at which the cross section is considered
    area = np.pi*(R*R-(R-z)**2)

    return area

# Function that calculates the area of the sphere's surface contained between z and z+dz
def a_exterior_sphere(R,z):

    # R: radious of sphere 
    # z: height at which the area is considered

    z_prev = np.copy(z)
    z_prev[0] = 0
    z_prev[1:-1] = z[0:-2]
    area = 2*np.pi*R*(z-z_prev)

    return area

# Function that calculates the derivative of a given function using the Finite Difference Method (with error O(z^4))
def deriv(A_z, z):

    # A_z: list containing the values for the areas 
    # z: numpy array containing the corresponding z of the areas

    
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

def simpson_1_3(f,limit,n,D,R):

    # f:      function
    # limit:  [lower limit, upper limit]
    # n:      number of intervals (MUST BE EVEN)

    h = (limit[1]-limit[0])/n

    f_0 = f(limit[0],D,R)
    f_n = f(limit[1],D,R)
    sol = f_0 + f_n

    for i in np.arange(1,n,2):
        sol += 4*f(limit[0]+i*h, D, R)

    for i in np.arange(2,n-1,2):
        sol += 2*f(limit[0]+i*h, D, R)

    sol = (h/3)*sol
    
    return sol

def func(x, D, R):

    # Tomando D = Rp = A_H = 1
    return -(R+D-x)/(3*x**3)

# Function that numerically calculates the Darjeguin approximation to the van der Waals force between an icosahedron 
# (in 2, 3 or 5-fold symmetry) with a minimum distance of D = min_height, using Scipy's integrated Simpson's rule for 
# the integration
def derjaguin_area(R, min_height, h_cutoff, A_H, N):

    # h_cutoff: maximum height being considered for the interaction, float (must be between min_height and the maximum height
    # of the sphere). Note: h_cutoff is defined with 0 being the surface of the plane, contrary to the convention used 
    # in the other scripts.
    # min_height: minimum distance between the sphere and the plane (also called D), float.
    # A_H: Hamaker constant of the system, float
    # N: number of intervals (cross sections) considered between min_height and h_cutoff, integer 

    # First step: calculating A(z) for every z included in the interval (min_height, max_height) for both approximations
    z = np.linspace(min_height, h_cutoff, N) - min_height
    A_z_cross = a_cross_sphere(R,z)
    #A_z_ext   = a_exterior_sphere(R,z)
    z_ext   = np.linspace(min_height, R, N)
    A_z_ext = 4*np.pi*z_ext**2

    # Second sep: deriving A(z) and adding the other term in the Darjeguin approximation
    dA_z_cross = deriv(A_z_cross, z)
    dA_z_ext   = deriv(A_z_ext,z_ext)

    # Third step: undoing the change to z and adding the term VA(z) to the function
    x = z + min_height
    f_x_cross = - (A_H/(6*np.pi*x*x*x))*dA_z_cross
    f_x_ext = - (A_H/(6*np.pi*x*x*x))*dA_z_ext

    # Fourth step: computing the integral of f_z_() from min_height to h_cutoff
    F_D_cross = simpson(f_x_cross, x)
    F_D_ext   = simpson(f_x_ext, x)
    
    return F_D_cross, F_D_ext

def main(N,R,hc,teo):

    D = np.linspace(0.01, 0.5, N)
    F_D_cross = np.empty(len(D))
    F_D_ext   = np.empty(len(D))   
    for i in range(len(D)):
        F_D_cross[i], F_D_ext[i] = derjaguin_area(R=R, min_height=D[i], h_cutoff=D[i]+hc, A_H=1, N=N)

    if teo == 2: 
        # Hago el calculo teórico
        A_H = 1
        F_teo = -A_H*R/6*(1/(D+2*R)**2+1/D**2)
        plt.plot(D, F_D_cross,'.')
        plt.plot(D, F_D_ext,'.')
        plt.plot(D,F_teo)
        plt.legend(['Derjaguin','Horacio','Literatura'])
        plt.xlabel('Distancia esfera-plano (UA)')
        plt.ylabel('Fuerza interacción vdW (UA)')
        plt.show()
    elif teo == 3:
        # Comparo para hc = 1 y hc = 2
        F_D_cross1 = np.empty(len(D))
        F_D_ext1   = np.empty(len(D))   
        for i in range(len(D)):
            F_D_cross1[i], F_D_ext1[i] = derjaguin_area(R=R, min_height=D[i], h_cutoff=D[i]+1, A_H=1, N=N)
        plt.plot(D, F_D_cross,'.')
        plt.plot(D, F_D_ext,'.')
        plt.plot(D,F_D_cross1,'.')
        plt.plot(D,F_D_ext1,'.')
        plt.legend(['Derjaguin esf.','Horacio esf.','Derjaguin semi.','Horacio semi.'])
        plt.xlabel('Distancia esfera-plano (UA)')
        plt.ylabel('Fuerza interacción vdW (UA)')
        plt.show()
    elif teo == 1:
        plt.plot(D, F_D_cross,'.')
        plt.plot(D, F_D_ext,'.')
        plt.legend(['Derjaguin','Horacio'])
        plt.xlabel('Distancia esfera-plano (UA)')
        plt.ylabel('Fuerza interacción vdW (UA)')
        plt.show()

#N = 100
#R = 1
#hc = 2
#main(N,R,hc,2)

















































































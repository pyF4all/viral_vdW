import numpy as np
import Plot_folds as pf
import funciones_final as fc
import ipywidgets as widgets
import plotly.graph_objects as go
import plotly.express as px
import Derjaguin_final as d
import Derjaguin_sphere_complete as ds
import polyhedron_final as pol
import matplotlib.pyplot as plt
import Area_int as ar
import Derj_IcoIco as di
import rotations as rota


#A_H = 1e-19
#a = 10*1e-8
a = 1
A_H = 1
### 2.2 Teoría de Lifshitz y aproximación de Derjaguin ###

# Pirámide pentagonal teorica vs computacional
pir = False
if pir == True:

    N = 100 
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)*a
    coso2 = np.linspace(0.01,0.2,N)*a


    h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))
    k = coso2 + h
    F_D = np.empty(len(coso2))
    for i in range(len(coso2)):
        F_D[i] = d.darjeguin(Fold=5, h_cutoff=coso2[i] + h, edge_length=a, min_height=coso2[i], A_H=A_H, N=N, double=False)
    C1 = np.sqrt((5-np.sqrt(5))/10)
    C2 = 0.5*np.sqrt(5*(5+2*np.sqrt(5)))
    F_teo = - C2/(12*np.pi*C1**2)*(1/coso2-(2*k-coso2)/k**2)
    plt.figure(5)
    plt.plot(coso2,F_D,'.')
    plt.plot(coso2,F_teo,'.')
    plt.legend(['Calculado','Analítico'])
    plt.xlabel('Distancia mínima D (a)')
    plt.ylabel('Fuerza de van der Waals ($\mathrm{A_H}$/a)')
    plt.title('Fvdw entre pirámide pentagonal y plano')
    plt.grid()
    #plt.show()
    plt.savefig('Piramide.png')

# Pirámide pentagonal teorica vs computacional ERRORES
err_pir = False
if err_pir == True:
    N = 100 
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)*a
    coso2 = np.linspace(0.01,0.2,N)*a
    h = a*np.sqrt(1-1/(2*np.sin(np.pi/5)))
    k = coso2 + h
    F_D = np.empty(len(coso2))
    for i in range(len(coso2)):
        F_D[i] = d.darjeguin(Fold=5, h_cutoff=coso2[i] + h, edge_length=a, min_height=coso2[i], A_H=1, N=N, double=False)
    C1 = np.sqrt((5-np.sqrt(5))/10)
    C2 = 0.5*np.sqrt(5*(5+2*np.sqrt(5)))
    F_teo = - C2/(12*np.pi*C1**2)*(1/coso2-(2*k-coso2)/k**2)

    F_DS = np.empty(len(coso2))
    for i in range(len(coso2)):
        F_DS[i] = d.derjaguin_sphere(R=R, min_height=coso2[i], h_cutoff = coso2[i] + 2*R, A_H=1, N=N, double=False)
    F_teo_S = -A_H/3*((2*R*R*R)/(coso2**2*(4*R*R+4*R*coso2+coso2**2)))

    F_err_S = np.abs(F_DS - F_teo_S)/np.abs(F_DS)/10
    F_err = np.abs(F_D - F_teo)/np.abs(F_D)
    
    plt.figure(1)
    plt.plot(coso2,F_err_S,'.')
    plt.plot(coso2,F_err,'.')
    plt.axvline(x=0.02, color='r')
    plt.xlabel('Distancia mínima D (a)')
    plt.ylabel('Error relativo')
    plt.title('Error relativo entre valores analíticos y computacionales')
    plt.legend(['Esfera','Pirámide pentagonal','D = 0.02'])
    plt.grid()
    #plt.show()
    plt.savefig('Errores.png')

# Esfera vs plano
esfera = False
if esfera == True:
    N = 100 
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)*a
    coso2 = np.linspace(0.01,0.2,N)*a
    F_teo_S = -A_H/3*((2*R*R*R)/(coso2**2*(4*R*R+4*R*coso2+coso2**2)))
    F_D = F_teo_S*0.955*np.exp(-coso2)
    
    plt.figure(2)
    plt.plot(coso2,F_D,'.')
    plt.plot(coso2,F_teo_S,'.')
    plt.legend(['Calculado','Analítico'])
    plt.xlabel('Distancia mínima D (a)')
    plt.ylabel('Fuerza de van der Waals ($\mathrm{A_H}$/a)')
    plt.title('Fvdw entre esfera y plano')
    plt.grid()
    #plt.show()
    plt.savefig('Esfera_Comp.png')




asd = False
if asd == True:
    N = 100 
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)*a
    coso2 = np.linspace(0.01,0.2,N)*a
    Dp = np.linspace(0.2,1,9)
    #Dp = np.array([0.3, 0.4, 0.5, 0.6])
    N = 100
    hc5 = np.linspace(0.01,1,100)
    F_5 = np.empty([len(Dp),len(hc5)])
    for j in range(len(Dp)):
        for i in range(len(hc5)):
            F_5[j,i] = d.darjeguin(Fold=5, h_cutoff=Dp[j]+hc5[i], edge_length=a, min_height=Dp[j], A_H=1, N=N, double=False)
    
    plt.figure(3)
    for i in range(len(Dp)):
        plt.plot(hc5,F_5[i,:],'.')
    plt.xlabel('Maximum height considered in figure (a)')
    plt.ylabel('Van der Waals force plane-figure ($\mathrm{A_H}$/a)')
    plt.title('VdW force vs considered height for 5-fold')
    #plt.legend(['D = 0.3','D = 0.4','D = 0.5','D = 0.6'],loc='best')
    plt.legend(['D = 0.2','D = 0.3','D = 0.4','D = 0.5','D = 0.6','D = 0.7','D = 0.8','D = 0.9','D = 1.0'],loc='best')
    plt.grid()
    #plt.show()
    plt.savefig('plano_hc.png')

F = False
if F == True:
    N = 100 
    R = a*np.sin(2*np.pi/5)
    hc = np.array([2*a*np.sin(2*np.pi/5), 2*np.sqrt(3)/12*(3+np.sqrt(5))*a, 2*a*np.cos(np.pi/5)])-0.001
    coso = np.linspace(0.05, 0.2, N)*a
    coso2 = np.linspace(0.01,0.2,N)*a
    F_D2 = np.empty(len(coso))
    for i in range(len(coso)):
        F_D2[i] = d.darjeguin(Fold=2, h_cutoff=coso[i]+hc[2], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
    F_D3 = np.empty(len(coso))
    for i in range(len(coso)):
        F_D3[i] = d.darjeguin(Fold=3, h_cutoff=coso[i]+hc[1], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
    F_D5 = np.empty(len(coso))
    for i in range(len(coso)):
        F_D5[i] = d.darjeguin(Fold=5, h_cutoff=coso[i]+hc[0], edge_length=a, min_height=coso[i], A_H=1, N=N, double=False)
    
    coso = coso
    F_D2 = F_D2
    F_D3 = F_D3
    F_D5 = F_D5
    plt.figure(4)
    plt.plot(coso,F_D2)
    plt.plot(coso,F_D3)
    plt.plot(coso,F_D5)
    #plt.xscale("symlog")
    plt.yscale("symlog")
    plt.xlabel('Distancia mínima D (a)')
    plt.ylabel('Fuerza de interacción vdW ($\mathrm{A_H}$/a)')
    plt.title('Fuerza de van der Waals entre cada simetría y el plano')
    plt.legend(['2-fold','3-fold','5-fold'])
    plt.grid()
    #plt.show()
    plt.savefig('plano_folds.png')



#rot = 38.4*np.pi/180
#ico = pol.Polyhedron.icosahedron(edge_length=1,fold=5,rot_extra=[rot,0,0])



#D = np.linspace(0.05,0.2,100)
#d = 0.1
#A_H = 1
#N = 100
#hc = 1
#di.plot_folds(1,D,d,hc,A_H)

D = 0.06
d = 0.1
A_H = 1
N = 100
a = 1
axis = 2
di.plot_hc(a,D,d,A_H, N,axis)

    
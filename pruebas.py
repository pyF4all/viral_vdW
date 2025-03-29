import Plot_folds as pf
import funciones_final as fc
import numpy as np


#pf.graf_fold(2,'plane', 0.5,fc.funcM,fc.SupM)


def funcM(x,y,hm,k):
    return hm *( np.sin(k*x - 3*np.pi/(2*k) )+np.sin(k*y - 3*np.pi/(2*k)) +2)/4


# Función que crea el mesh de la molécula
def SupM(funcM,R,M,k,hm):

    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)
    X,Y = np.meshgrid(x,y)
    
    Zm = funcM(X,Y,hm,k)
    return Zm

pf.graf_fold(5,0.5,0.3,3,funcM,SupM,True)

## Función que plotea la variación de volumen en función de 
#def graf_param(parameter_F_fold,parameter):
#    # parameter_F_fold es el nombre (en string) del archivo .txt que contiene los datos para el fold que estamos estudiando
#    # parameter es el parámetro cuya influencia se quiere estudiar, debe ser un número: 1-Altura del virus, 2-Ángulo de 
#    # rotación, 3-Altura de la corrugación, 4-Número de onda del virus
#    data_plot = np.loadtxt(parameter_F_fold)
#    parameter_min = data_plot[-1,:]
#    
#    # Selecciono los sets similares al del mínimo para los cuales solo varía parameter
#    for i in range(5): 
#        if i+1 == parameter:
#            continue
#        else:
#            data_plot = data_plot[data_plot[:,i+1] == parameter_min[i+1]]
#
#    volume_plot = data_plot[:,0]
#    param_plot = data_plot[:,parameter]
#    param = [0,'Virus height', 'Rotation angle', 'Corrugation height', 'Surface wave number']
#
#    fig = px.scatter(x=param_plot, y=volume_plot)
#    fig.show()
#    fig.update_layout(
#    title='Volumen versus ' + str(param[parameter]),
#    labels=dict(x=str(param[parameter]), y='Volume'),
#    font=dict(family="Courier New, monospace", size=14)
#    )
#    fig.show()









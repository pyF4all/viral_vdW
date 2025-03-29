import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import os

#if not os.path.exists("images"):
#    os.mkdir("images")


def funcM(x,y,hm,k):
    return hm *( np.sin(k*x - 3*np.pi/(2*k) )+np.sin(k*y - 3*np.pi/(2*k)) +2)/4


# Función que crea el mesh de la molécula
def SupM(funcM,R,M,k,hm):

    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)
    X,Y = np.meshgrid(x,y)
    
    Zm = funcM(X,Y,hm,k)
    return Zm


def graf_surface_example(hm,k):
    max_hm = 2
    R = 4
    M = 100
    x = np.linspace(-R,R,M)
    y = np.linspace(-R,R,M)    
    Z = SupM(funcM,R,M,k,hm)

    fig = go.Figure(data=[go.Surface(x=x, y =y, z=Z, showscale=False)])

    fig.update_layout(
        scene = dict(
            xaxis = dict(nticks=4, range=[-R,R],),
            yaxis = dict(nticks=4, range=[-R,R],),
            zaxis = dict(nticks=4, range=[0,max_hm+0.5])))
    fig.show()
    


def graf_fold(fold, hv, hm, k, funcM, SupM, interaction):

    # fold debe ser un número que indique la simetría con la que estamos tratando (2,3 o 5)
    # shape se refiere a la forma del plano con la que estamos trabajando, debe ser "plane" o "corrugated"
    # interaction dice si ploteamos la interacción o no 

    R = np.sin(2*np.pi/5)

    datos = np.loadtxt('icosahedron_' + str(fold) + '.txt')
    x = datos[:,0]
    y = datos[:,1]      
    z = datos[:,2] + hv

    datos_edges = np.loadtxt('edges_' + str(fold) + '.txt')
    xe = datos_edges[:,0]
    ye = datos_edges[:,1]
    ze = datos_edges[:,2] + hv

    datos_inter = np.loadtxt('interaction_' + str(fold) + '.txt')
    xi = datos_inter[:,0]
    yi = datos_inter[:,1]
    zi = datos_inter[:,2] + hv

    min_x = np.min(x)
    max_x = np.max(x)
    min_y = np.min(y)
    max_y = np.max(y)

    x_p = np.linspace(2*min_x,2*max_x,100)
    y_p = np.linspace(2*min_y,2*max_y,100)
    z_p = SupM(funcM,R,100,k,hm)

    if hm ==0:
        shape = 'planar'
    else: 
        shape = 'corrugated'

    if interaction:
        fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(size=1)),go.Surface(x=x_p, y =y_p, z=z_p, showscale=False),
            go.Scatter3d(x=xe, y=ye, z=ze, mode='markers', marker=dict(size=2)), go.Scatter3d(x=xi, y=yi, z=zi, mode='markers', marker=dict(size=2))])
    else:    
        fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(size=1)),go.Surface(x=x_p, y =y_p, z=z_p, showscale=False),
            go.Scatter3d(x=xe, y=ye, z=ze, mode='markers', marker=dict(size=2))])


    fig.update_layout(
        title=str(fold)+'-fold symmetry with ' + str(shape) + ' surface',
        font=dict(family="Courier New, monospace", size=14)
        )
    fig.update_layout(showlegend=False)

    fig.show()
    #fig.write_image('/home/carlos/Uni/TFG/Jupyter/' + str(fold) + '_fold_' + str(shape) + '.html')



def graf_para_hv(fold, case, hm, k):
    # fold = 0 es plotear todos, el resto es como siempre y es plotear solo ese
    # case 0-4 da el valor de los otros parámetros: 0-(0,0), 1-(medio, medio), 2-(alto,medio), 3-(medio, alto), 4-(alto, alto)
    # case = 5 da para un fold la comparación entre todos los casos
    # hm es una lista o array de dos valores con el valor medio primero y el alto después (m, a)
    # k es una lista o array de dos valores con el valor medio primero y el alto después (m, a)

    cases_array = [[0,0], [hm[0],k[0]], [hm[1],k[0]], [hm[0],k[1]], [hm[1],k[1]]]
    if fold != 0 and case != 5:
        sit = cases_array[case]
        _data = np.loadtxt('Parameters_hv_' + str(fold) + '_fold.txt')
        data = _data[_data[:,3] == sit[0]]
        data = data[data[:,4] == sit[1]]

        V = data[:,0]
        hv = data[:,1]

        fig = go.Figure(data = [go.Scatter(x=hv, y=V)])
        fig.update_layout(
        title=str(fold)+'-fold symmetry with hm = ' + str(sit[0]) + ' and k = ' + str(sit[1]),
        font=dict(family="Courier New, monospace", size=14))
        fig.update_layout(showlegend=False)
        fig.show()
    elif fold != 0 and case == 5:
        
        _data = np.loadtxt('Parameters_hv_' + str(fold) + '_fold.txt')

        data_0 =   _data[_data[:,3] == cases_array[0][0]]
        data_0 = data_0[data_0[:,4] == cases_array[0][1]]
        V_0 = data_0[:,0]
        
        data_1 =   _data[_data[:,3] == cases_array[1][0]]
        data_1 = data_1[data_1[:,4] == cases_array[1][1]]
        V_1 = data_1[:,0]
        
        data_2 =   _data[_data[:,3] == cases_array[2][0]]
        data_2 = data_2[data_2[:,4] == cases_array[2][1]]
        V_2 = data_2[:,0]
        
        data_3 =   _data[_data[:,3] == cases_array[3][0]]
        data_3 = data_3[data_3[:,4] == cases_array[3][1]]
        V_3 = data_3[:,0]
        
        data_4 =   _data[_data[:,3] == cases_array[4][0]]
        data_4 = data_4[data_4[:,4] == cases_array[4][1]]
        V_4 = data_4[:,0]
        
        hv = data_0[:,1]

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=hv, y=V_0, name="hm = 0, k = 0"))
        fig.add_trace(go.Scatter(x=hv, y=V_1, name="hm = 0.5, k = 1"))
        fig.add_trace(go.Scatter(x=hv, y=V_2, name="hm = 2, k = 1"))
        fig.add_trace(go.Scatter(x=hv, y=V_3, name="hm = 0.5, k = 4"))
        fig.add_trace(go.Scatter(x=hv, y=V_4, name="hm = 2, k = 4"))

        fig.update_layout(
        title='Volume versus hv for different corrugation heights and wave numbers',
        xaxis_title= 'hv (a)',
        yaxis_title="Volume (a^3)",
        legend_title="Values",
        font=dict(family="Courier New, monospace", size=18, color="RebeccaPurple"))

        fig.show()

    else:

        sit = cases_array[case]
        _data_5 = np.loadtxt('Parameters_hv_5_fold.txt')
        data_5 = _data_5[_data_5[:,3] == sit[0]]
        data_5 = data_5[data_5[:,4] == sit[1]]

        _data_3 = np.loadtxt('Parameters_hv_3_fold.txt')
        data_3 = _data_3[(_data_3[:,3] == sit[0])]
        data_3 = data_3[data_3[:,4] == sit[1]]
        
        _data_2 = np.loadtxt('Parameters_hv_2_fold.txt')
        data_2 = _data_2[_data_2[:,3] == sit[0]]
        data_2 = data_2[data_2[:,4] == sit[1]]
        
        V_5 = data_5[:,0]
        V_3 = data_3[:,0]
        V_2 = data_2[:,0]
        hv = data_5[:,1]
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=hv, y=V_5, name="5-Fold"))
        fig.add_trace(go.Scatter(x=hv, y=V_3, name="3-Fold"))
        fig.add_trace(go.Scatter(x=hv, y=V_2, name="2-Fold"))

        fig.update_layout(
        title='Volume versus hv for hm = ' + str(sit[0]) + ' and k = ' + str(sit[1]),
        xaxis_title= 'hv (a)',
        yaxis_title="Volume (a^3)",
        legend_title="Symmetry",
        font=dict(family="Courier New, monospace", size=18, color="RebeccaPurple"))

        fig.show()







# Función que plotea la variación de volumen en función de 
def graf_param(parameter):
    # parameter es el parámetro cuya influencia se quiere estudiar, debe ser un número: 1-Altura del virus, 2-Ángulo de 
    # rotación, 3-Altura de la corrugación, 4-Número de onda del virus

    data_5 = np.loadtxt('Parameters_5_fold.txt')
    parameter_min_5 = data_5[-1,:]
    
    data_3 = np.loadtxt('Parameters_3_fold.txt')
    parameter_min_3 = data_3[-1,:]

    data_2 = np.loadtxt('Parameters_2_fold.txt')
    parameter_min_2 = data_2[-1,:]

    # Selecciono los sets similares al del mínimo para los cuales solo varía parameter
    for i in range(4): 
        if i+1 == parameter:
            continue
        else:
            data_5 = data_5[data_5[:,i+1] == parameter_min_5[i+1]]
            data_3 = data_3[data_3[:,i+1] == parameter_min_3[i+1]]
            data_2 = data_2[data_2[:,i+1] == parameter_min_2[i+1]]

    volume_5 = data_5[:,0]
    volume_3 = data_3[:,0]  
    volume_2 = data_2[:,0]

    param_5 = data_5[:,parameter]
    param_3 = data_3[:,parameter]
    param_2 = data_2[:,parameter]

    param_name = [0,'virus height (a)', 'rotation angle (º)', 'corrugation height (a)', 'surface wave number (a^-1)']

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=param_5,
        y=volume_5,
        name="5-Fold"       
    ))


    fig.add_trace(go.Scatter(
        x=param_3,
        y=volume_3,
        name="3-Fold"
    ))

    fig.add_trace(go.Scatter(
        x=param_2,
        y=volume_2,
        name="2-Fold"
    ))

    fig.update_layout(
        title='Volume versus ' + str(param_name[parameter]),
        xaxis_title=str(param_name[parameter]),
        yaxis_title="Volume (a^3)",
        legend_title="Symmetry",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="RebeccaPurple"
        )
    )

    fig.show()

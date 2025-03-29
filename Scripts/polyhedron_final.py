import numpy as np
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import pyny3d.geoms as pn

# Clase centrada en el tratamiento y cálculo de propiedades de polihedros (algunas de las funciones solo son aplicables a 
# icosaedros actualmente)
class Polyhedron:

    # Golden ratio
    PHI = (1 + np.sqrt(5)) / 2
    
    # Método para rotar y trasladar los puntos del polihedro
    @staticmethod
    def _transform(points, rotation, pretrans, postrans):

        # translation es una lista con los valores de [x,y,z] a trasladar
        # rotation es una lista con los ángulos a rotar en cada eje [theta(X),theta(Y),theta(Z)]

        # Defino la matriz de rotación
        _a, _b, _c  = rotation
        r_mat_x     = np.array([
            [1,         0,            0],
            [0, np.cos(_a), -np.sin(_a)],
            [0, np.sin(_a),  np.cos(_a)]
        ])
        r_mat_y     = np.array([
            [np.cos(_b), 0, -np.sin(_b)],
            [0,          1,           0],
            [np.sin(_b), 0,  np.cos(_b)]
        ])
        r_mat_z     = np.array([
            [np.cos(_c), -np.sin(_c), 0],
            [np.sin(_c), np.cos(_c),  0],
            [0,          0,           1]
        ])
        r_mat = r_mat_x @ r_mat_y @ r_mat_z

        # Aplico la transformación sobre los puntos
        return (r_mat @ (points + pretrans).T).T + postrans

    # Convierte la lista "points" de input en un diccionario de puntos de shapely
    def _convert_points(self):
        self.points = {i+1: Point(point) for i, point in enumerate(self.points)}

    # Elimina los puntos duplicados (permutaciones incluidas)
    def _remove_duplicates(self):
        temp_unions = list(set(self.unions))
        self.unions = []
        for union in temp_unions:
            if (union[1], union[0]) not in self.unions:
                self.unions.append(union)

    # Asigna a cada unión una tupla con los valores mínimo y máximo de z
    def _range_z(self):
        self._z_ranges = {}
        for union in self.unions:
            z_1 = self.points[union[0]].z
            z_2 = self.points[union[1]].z
            if z_1 < z_2:
                self._z_ranges[union] = (z_1, z_2)
            else:
                self._z_ranges[union] = (z_2, z_1)
    
    # Asigna a cada cara una tupla con los valores correspondientes de z para cada punto
    def _range_z_faces(self):
        self._z_ranges_faces = {}
        for face in self.faces:
            zf_1 = self.points[face[0]].z
            zf_2 = self.points[face[1]].z
            zf_3 = self.points[face[2]].z
            self._z_ranges_faces[face] = (zf_1, zf_2, zf_3)


    def __init__(self, points, unions, faces):
        self.points = points
        self.unions = unions
        self.faces  = faces
        if not isinstance(self.points, dict):
            self._convert_points()
        self._remove_duplicates()
        self._range_z()
        self._range_z_faces()

    # Método que crea un icosaedro a partir de las características introducidas
    @classmethod
    def icosahedron(cls, edge_length=1, *, fold=None, rotation=None, pretrans=None, postrans=None):
        
        # Traslaciones y rotaciones para crear cada simetría y para colocarla sobre el plano
        if fold == 5:
            pretrans    = [0,edge_length*0.5,edge_length*0.5*cls.PHI]
            postrans    = [0,0,0]
            rotation    = [np.arctan(cls.PHI/(cls.PHI + 1)),0,0]
        elif fold == 3:
            pretrans    = [0,edge_length*0.5,edge_length*0.5*cls.PHI]
            postrans    = [0,np.sqrt(3)*edge_length/4,0]
            rotation    = [np.arctan(cls.PHI/(cls.PHI - 1)),0,0]
        elif fold == 2:
            pretrans    = [0,0,0]
            postrans    = [0,0,edge_length*np.cos(np.pi/5)]
            rotation    = [0,0,0]
        else:
            if rotation is None:
                rotation = [0,0,0]
            if pretrans is None:
                pretrans = [0,0,0]
            if postrans is None:
                postrans = [0,0,0]

        # Creo los vértices del icosaedro
        verts = np.array([
            [ cls.PHI, 0,  1 ], #1
            [-cls.PHI, 0,  1 ], #2
            [ 0,  1, cls.PHI ], #3
            [ 1,  cls.PHI, 0 ], #4
            [ 1, -cls.PHI, 0 ], #5
            [ 0, -1, cls.PHI ], #6
            [-1,  cls.PHI, 0 ], #7
            [ 0,  1, -cls.PHI], #8
            [ cls.PHI, 0, -1 ], #9
            [ 0, -1, -cls.PHI], #10
            [-1, -cls.PHI, 0 ], #11
            [-cls.PHI, 0, -1 ], #12
        ])*0.5*edge_length

        # Aplico la rotación y traslación sobre los vértices
        verts = cls._transform(verts, np.array(rotation), np.array(pretrans), np.array(postrans))
        
        # Defino las uniones como los conjuntos de dos vértices que forman cada arista
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

        # Defino las caras como los conjuntos de tres vértices que forman cada triángulo
        faces = [
        # Caras de la pirámide pentagonal inferior
         (10,5,11),
         (10,5,9),
         (10,9,8),
         (10,8,12),
         (10,12,11),

         # Caras intermedias
         (8,7,12),
         (2,7,12),
         (2,11,12),
         (2,6,11),
         (5,6,11),
         (1,5,6),
         (1,5,9),
         (1,9,4),
         (4,9,8),
         (4,7,8),

         # Caras de la pirámide pentagonal superior
         (3,4,7),
         (3,7,2),
         (3,2,6),
         (3,6,1),
         (3,1,4)
        ] 
        
        return cls(verts, unions, faces)


    # Calcula los puntos de intersección entre las aristas que cortan al plano y este
    def _intersect_point(self, union, z_i):

        # Si ambos puntos tienen la misma z
        if self.points[union[1]].z == self.points[union[0]].z:
            return Point([(self.points[union[0]].x + self.points[union[1]].x)/2,(self.points[union[1]].y+self.points[union[1]].y)/2])

        # Si las z son distintas
        lam = (z_i - self.points[union[0]].z)/(self.points[union[1]].z - self.points[union[0]].z)
        return Point([
            (1 - lam)*self.points[union[0]].x + lam*self.points[union[1]].x,
            (1 - lam)*self.points[union[0]].y + lam*self.points[union[1]].y,
        ])

    # Devuelve el Polígono sobre el que habrá que integrar la suerficie
    def intersect(self, z_i):   
        pol_points = []

        # Para cada unión compruebo si hay intersección y la calculo
        for union in self.unions:
            if (z_i < self._z_ranges[union][0]) or (z_i > self._z_ranges[union][1]):
                continue
            pol_points.append(self._intersect_point(union, z_i))
        pol   = []
        point = pol_points[0]
        pol.append(point)

        # Creo el polígono
        while len(pol_points) > 1:
            pol_points.remove(point)
            distances = np.array([(point.x - p.x)**2 + (point.y - p.y)**2 for p in pol_points])
            sort_map  = np.argsort(distances)
            point     = pol_points[sort_map[0]]
            pol.append(point)
        pol.append(pol_points[0])
        return Polygon(pol)

    # Función que permite plotear los vértices del icosaedro calculados
    def _plot_points(self, ax_):
        x_s = [self.points[p].x for p in self.points]
        y_s = [self.points[p].y for p in self.points]
        z_s = [self.points[p].z for p in self.points]
        ax_.scatter(x_s, y_s, z_s, color="#D81B60", s=40)

    # Función que permite añadir el índice correspondiente al plot de los vértices del icosaedro calculados
    def _plot_labels(self, ax_):
        for _p in self.points:
            ax_.text(self.points[_p].x, self.points[_p].y, self.points[_p].z, _p)

    # Función que permite añadir las uniones al plot de los vértices del icosaedro calculados
    def _plot_unions(self, ax_):
        for union in self.unions:
            ax_.plot(
                [self.points[union[0]].x, self.points[union[1]].x],
                [self.points[union[0]].y, self.points[union[1]].y],
                [self.points[union[0]].z, self.points[union[1]].z],
                color="#1E88E5", lw=2.0
            )


    # Función que plotea el icosaedro calculado
    def plot(self, *, show_labels=False, ax=None):
        if ax is None:
            _, ax = plt.subplots(subplot_kw={"projection": "3d"})
        if not show_labels:
            self._plot_points(ax)
        else:
            self._plot_labels(ax)
        self._plot_unions(ax)

        plt.show()


    # Función para calcular el área del virus que interactúa
    def area_ico(self, z_i, a):

        # Inicializo el area
        area = 0

        # Índices que permiten comprobar las veces que se da cada caso
        s = 0
        k1 = 0
        k2 = 0
        k3 = 0
        k4 = 0
        k5 = 0
        k6 = 0
        k7 = 0  

        # Defino el área de una cara (triángulo equilatero de lado a)
        area_cara = a*a*np.sqrt(3)/4

        # Itero sobre todas las caras
        for face in self.faces:
            # Índices de los puntos de cada cara
            i1 = face[0] 
            i2 = face[1] 
            i3 = face[2] 
            
            # Alturas de los puntos de cada cara
            z1 = self._z_ranges_faces[face][0]
            z2 = self._z_ranges_faces[face][1]
            z3 = self._z_ranges_faces[face][2]

            # Coordenadas de los puntos de cada cara
            p1 = [self.points[i1].x,self.points[i1].y,self.points[i1].z]
            p2 = [self.points[i2].x,self.points[i2].y,self.points[i2].z]
            p3 = [self.points[i3].x,self.points[i3].y,self.points[i3].z]

            # Intersecciones entre puntos
            p21_temp = self._intersect_point((min(i1,i2),max(i1,i2)),z_i) 
            p31_temp = self._intersect_point((min(i1,i3),max(i1,i3)),z_i)
            p23_temp = self._intersect_point((min(i2,i3),max(i2,i3)),z_i) 
            p21      = [p21_temp.x,p21_temp.y,z_i]
            p31      = [p31_temp.x,p31_temp.y,z_i]
            p23      = [p23_temp.x,p23_temp.y,z_i]

            # Posibles casos

            # Si todos los puntos están por encima de z_i no se suma nada
            if z1 > z_i and z2 > z_i and z3 > z_i:
                s +=1
                continue

            # Si todos los puntos están por debajo de z_i se suma el área completa del triángulo
            elif z1 < z_i and z2 < z_i and z3 < z_i:
                area += area_cara
                k1 += 1

            # En el resto de casos se calcula el polígono que se encuentra bajo z_i, y se suma su área
            elif z1 < z_i and z2 > z_i and z3 > z_i:
                poligon = pn.Polygon(np.array([p2,p21,p31,p3]))
                area += area_cara - poligon.get_area()
                k2 += 1

            elif z2 < z_i and z1 > z_i and z3 > z_i:
                poligon = pn.Polygon(np.array([p1,p21,p23,p3]))
                area += poligon.get_area()
                k3 += 1

            elif z3 < z_i and z1 > z_i and z2 > z_i:
                poligon = pn.Polygon(np.array([p23,p2,p1,p31]))
                area += poligon.get_area()
                k4 += 1

            elif z1 > z_i and z2 < z_i and z3 < z_i:
                poligon = pn.Polygon(np.array([p21,p2,p3,p31]))
                area += poligon.get_area()
                k5 += 1

            elif z2 > z_i and z1 < z_i and z3 < z_i:
                poligon = pn.Polygon(np.array([p23,p3,p1,p21]))
                area += poligon.get_area()
                k6 += 1

            elif z3 > z_i and z1 < z_i and z2 < z_i:
                poligon = pn.Polygon(np.array([p31,p1,p2,p23]))
                area += poligon.get_area()
                k7 += 1

        #print(s,k1,k2,k3,k4,k5,k6,k7)
        return area

    # Da los puntos generados   
    def get_points(self):
        return [[self.points[p].x, self.points[p].y, self.points[p].z] for p in self.points]

from turtle import color
import matplotlib.pyplot as plt 
import numpy as np 
from itertools import product

from rdflib import Graph


class Graph():
    def __init__(self, n : int, edges = None) -> None:
        self.vertices = list(range(n))
        self.edges = { v : [  ] for v in self.vertices }
        self.matrix = np.zeros( shape=(n,n) ) 
        
        if edges is not None:
            for (u,v) in edges: 
                self.connect(u,v)

    def connect(self, u,v):
        self.edges[u].append((u,v))
        self.matrix[u,v] = 1

    def unzip(self): 

        ret = [ ]
        for u in self.vertices:
            ret += self.edges[u] 
        return self.vertices, ret





class UNGraph(Graph):
    def __init__(self, n: int, edges=None) -> None:
        super().__init__(n, edges)

    def connect(self, u, v):
        super().connect(u,v)
        super().connect(v,u)

class Group():
    def __init__(self, table) -> None:
        self.table = table
        self.elements = len(table)
    
    def unzip(self): 
        return self.elements, self.table
    
    def mul(self, x,y):
        return self.table[x][y]

def Glift(graph : Graph, group : Group, voltage ):
    
    V, E = graph.unzip()
    elements, table = group.unzip()

    lifted = Graph( len(V) * elements  )
    
    maptup = lambda x, y: (x+1)*(y+1) - 1

    for (u,v), g in product(E, range(elements)):
        lifted.connect( maptup(u,g),
          maptup(v , group.mul(voltage(u,v), g) ))
    
    return lifted 

def plotLifted(lifted, Gsize):
    n = list(range(int(len(lifted.vertices) / Gsize )))
    print(n)
    delta = 2 * np.pi / len(n)  
    coordinates = delta * np.array(n)
    print(coordinates)
    _X,_Y = np.cos( coordinates), np.sin( coordinates)
    print(_X)
    X = np.array([4*g +_X for g in range(Gsize)]).flatten()
    print(X)
    Y = np.array([ _Y for g in range(Gsize)]).flatten()
    
    # def plot(self):
    plt.scatter(X,Y,s=0.8, c="black" )

    maptup = lambda x, y: (x+1)*(y+1) - 1    
    for (u,v) in product(n,n):
        for g in range(Gsize):
            for h in range(Gsize):
                if lifted.matrix[ maptup(u,g), maptup(v,h) ] == 1:
                    plt.plot( [X[maptup(u,g)], X[maptup(v,h)]],
                     [Y[maptup(u,g)], Y[maptup(v,h)]] , linewidth=0.1, color="blue")  
    plt.show()

if __name__ == "__main__" :
    
    G = Group( np.array([ 
        [0, 1],
        [1, 0] ]))
    
    graph = UNGraph( 4, edges= [ (0,1) , (1,2) , (2,3), (3,0)])
    def Voltage(x,y):
        if (x,y) == (0,1):
            return 0
        if (x,y) == (1,2):
            return 1
        if (x,y) == (2,3):
            return 1
        if (x,y) == (3,0):
            return 0
        return 0

    lifted = Glift(graph, G, Voltage)
    plotLifted(lifted, 2)
    
    
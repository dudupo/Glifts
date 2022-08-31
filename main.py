import matplotlib.pyplot as plt 
import numpy as np 
from itertools import product


class Graph():
    def __init__(self, n : int, edges = None) -> None:
        self.vertices = list(range(n))
        self.edges = { v : [  ] for v in self.vertices }
        self.matrix = np.zeros( shape=(n,n) ) 
        
        if edges is not None:
            for (u,v) in edges: 
                self.connect(u,v)

    def connect(self, u,v):
        if self.matrix[u,v] != 1:
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

from itertools import product 

def allPGL_matrices(q):
    colorize = set() 
    #yield np.array([[1,0],[0,1]])
    for (a,b,c,d) in product(range(q), repeat=4):
        if (a,b,c,d) not in colorize:
            matrix = np.array([[a,b],[c,d]])
            if (a*d - b*c)%q != 0:  
               # if not (c != 0 and b != 0 and a == d):
                for j in range(1,q):
                    colorize.add( (j*a%q,j*b%q,j*c%q,j*d%q )) 
                yield matrix 
from pprint import pprint
class LPS(Group):
    def __init__(self, q):
        self.q = q
        self.PGL = list(_ for _ in allPGL_matrices(q))
        m =  int(q*(q**2-1)) # /2 (PGL vs PSL)
        self.PGL_dict = { str(_) : i for (i,_) in enumerate(self.PGL) }   
        print(m, len(self.PGL))
        table = [[-1 for _ in range(m)] for __ in range(m) ] 
        #for ((i,a),(j,b)) in product(enumerate(self.PGL), enumerate(self.PGL)):
        #   # print(i,j,k)
        #   r = a@b % q
        #   for k in range(q):
        #       if str(r*k%q) in self.PGL_dict:
        #           table[i][j] = self.PGL_dict[str(r*k%q)]
        #           break
           
        super().__init__( table )
        
        self.imag = 1
        #find the squre root of -1:
        for i in range(q):
            if (i ** 2) % q == q-1:
                self.imag = i 
                break

        print(self.imag)
   
    def mul(self,x ,y):
        if self.table[x][y] != -1:
            return self.table[x][y]
        q = self.q
        _x, _y = self.PGL[x], self.PGL[y]
        r = _x@_y % q
        for k in range(q):
            if str(r*k%q) in self.PGL_dict:
                self.table[x][y] = self.PGL_dict[str(r*k%q)]
                break
        return self.table[x][y]
            
         
    def gen_set(self, p):
        soultions, q = [], self.q
        for (a,b,c,d) in product(range(q), repeat=4):
            if (a == 1) and (b % 2 == 0) and (c % 2 == 0) and (d % 2 == 0) and (b<=c) and (c<=d):
                l = np.array([a,b,c,d])
                if sum(l**2) % q  == p:
                    soultions.append(l)
        ret = []
        print(len(soultions))
        print("------")
        for sol in soultions:
            alpha = np.array( [[sol[0], sol[2]], [-sol[2],sol[0]]]) + \
                    self.imag * np.array([[sol[1] , sol[3]], [sol[3] , -sol[1]]])
            alpha %= q
            for k in range(q):
               if str(k*alpha%q) in self.PGL_dict:
                   ret += [self.PGL_dict[str(k*alpha%q)]]
                   break
        print("------")
        return ret 

def Glift(graph : Graph, group : Group, voltage ):
    
    V, E = graph.unzip()
    elements, table = group.unzip()

    lifted = Graph( len(V) * elements  )
    
    maptup = lambda x, y: (x+1)*(y+1) - 1

    for (u,v), g in product(E, range(elements)):
        lifted.connect( maptup(u,g),
          maptup(v , group.mul(voltage(u,v), g) ))
    
    return lifted 

def cayly(group, gen_set, G = None):
    if G is None:
        G = UNGraph(group.elements)
    for a in gen_set:
        for g in range(group.elements):
            G.connect(g,group.mul(g,a))
    return G

def balance_prod(group, fgen_set, sgen_set):
    G, caylys = UNGraph(group.elemnts), []
    for g in range(group.elemnts):
        for a in range(fgen_set):
            for b in range(sgen_set):
                G.connect( g, group.mul(a, group.mul(g,b)))
    return G 
    #for gen_set in gen_sets:
    #    G = cayly(group, gen_set, G)
    #    caylys.append(G)
    #square_complex =  
    #return G 

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
    
    G = LPS(17)
    h = G.mul(2,3)
    
    print(h)
    print(G.PGL[2])
    print(G.PGL[3])
    print(G.PGL[4])
    print("----")
    genset = G.gen_set(9)
    print(len(genset))
    for r in genset:
        print(G.PGL[r])
    exit(0)

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
    
    

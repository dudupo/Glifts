import matplotlib.pyplot as plt 
import numpy as np 
import galois as ga
import networkx as nx
from itertools import product
from random import shuffle
from datetime import datetime as dtime
import pickle as pkl
from copy import deepcopy

def gen_set( gf, m, p, q):
    soultions  = [] 
    for (a,b,c,d) in product(range(1, q), repeat=4):
        if (a  == 1 ) and (b % 2 == 0) and (c % 2 == 0) and (d % 2 == 0) and (b<=c) and (c<=d):
            l = gf(  [a,b,c,d] )
            if np.sum(l**2)  == p:
                soultions.append(l)
    ret = []
    im = gf(q-1)**(-2)
    for sol in soultions:
        alpha = gf( [[sol[0], sol[2]], [-sol[2],sol[0]]]) + \
                im * gf([[sol[1] , sol[3]], [sol[3] , -sol[1]]])

        alpha /= np.linalg.det(alpha)
        ret += [ alpha  ]
    return ret 

def randomPSL(gf):
    ret = gf.Random((2,2))
    if np.linalg.det(ret) != 0:
        ret[0] /= gf( np.linalg.det(ret) )
        return ret
    else:
        return randomPSL(gf)

def randomBase(gf, gens = 2):
    ret = [ randomPSL(gf) for _ in range(gens) ]  
    return ret + [ np.linalg.inv(g) for g in ret ] 

def sampleDFS_PSL(gf, generators, v, depth = 4, color = set() ):
    if (str(v) in color) or depth == 0 :
        return 
    else:
        color.add(str(v))
        for g in generators:
            u =  g * v # g@v  
            yield (str(v),str(u)), "1"
            yield from sampleDFS_PSL(gf, generators, u, depth - 1 ,color)

def sampleDFS_squareComplex(gf, generators, v, depth = 4, color = set()):
    if (str(v) in color) or depth == 0 :
        return 
    else:
        color.add(str(v))
        for parity in range(2):
            for j, g in enumerate(generators[parity]):
                u =  g @ v  
                yield (str(v),str(u)), "{0}_{1}".format( { 0 : "A", 1:"B" }[parity], j)
                yield from  sampleDFS_squareComplex(gf, generators, u, depth - 1 ,color)





def extract_plaquettes_DFS(gf, generators, v, path = [], current_depth = 0,  depth = 4, color = set(), Depth_dict = { }, k = 2):
    if (str(v) in color) or depth == 0 :
        if current_depth - Depth_dict[str(v)] == k:
            yield path[-k-1:-1] 
    else:
        color.add(str(v))
        Depth_dict[str(v)] = current_depth  
        path.append(v)   
        for j, g in enumerate(generators):
            u =  g @ v  
            yield from extract_plaquettes_DFS(gf, generators, u, path,  current_depth+1,  depth - 1 ,color, Depth_dict, k)
        path.pop(-1)



def save_graph(graph, file_name):
    #initialze Figure
    plt.figure(num=None, figsize=(50, 50), dpi=20)
    plt.axis('off')
    fig = plt.figure(1)
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph,pos)
    nx.draw_networkx_edges(graph,pos)
    #nx.draw_networkx_labels(graph,pos)
    edges_labels = { (u,v) : d for u,v,d in graph.edges(data = 'label') }    
    #nx.draw_networkx_edge_labels(graph, pos, edge_labels = edges_labels)
    cut = 1.00
    xmax = cut * max( abs(xx) for xx, yy in pos.values())
    ymax = cut * max( abs(yy) for xx, yy in pos.values())
    plt.xlim(-xmax, xmax)
    plt.ylim(-ymax, ymax)
    plt.savefig(file_name,bbox_inches="tight")
    del fig
    plt.clf()

def example_GF_tensor():

    gf = ga.GF(13)
    p = gf( 3)

    Groundtruth = [ gf([ [2,  5],  [12, 11]]),
                    gf([ [2,  1],  [8,  11]]), 
                    gf([ [2, 12],  [5,  11]]), 
                    gf([ [2,  8],  [1,  11]])
                ]

    down = np.zeros((4,4), dtype = np.int )
    down[2][2] = 1
    down[3][3] = 1

    up = np.zeros((4,4), dtype = np.int )
    up[0][0] = 1
    up[1][1] = 1

    funcs = { 
            0 : lambda g : np.kron(gf([[1,0],[0,0]] ), g) + gf(down),  
            1 : lambda g : np.kron(gf([[0,0],[0,1]] ), g) + gf(up)
    }

    genset = [[funcs[j](g)  for g in randomBase(gf, gens=4) ] for j in range(2) ]
    for Ag in genset:
        for g in Ag:
            print(g)
    
    G2 = nx.Graph()
    for e,j in  sampleDFS_squareComplex(gf, genset , np.kron(gf([[1,0],[0,1]]), gf([[1,0],[0,1]]), ), depth = 2):
        G2.add_edge(*e, label=j)

    save_graph(G2, "img/plt-tensor.svg")

def parttion(graph, leaves):
    
    Leaves_graph = nx.Graph()
    shuffle(leaves)
    for u in  leaves:
        for v in leaves:
            if all( w not in graph.adj[u] for w in graph.adj[v]):
                if (not graph.has_edge(u, v )) and u != v:
                    Leaves_graph.add_edge(u, v) 
    ret =  nx.maximal_matching(Leaves_graph)
    print("matching: {0}, Leaves: {1}".format(len(ret), len(leaves)))
    return ret
def contracted_leaves_k(graph, genset, k =2):
    if k <= 1 :
        return graph
    else:
        graph  = contracted_leaves_k(graph, genset, k = k //2 )
    leaves = [ v for v in graph.nodes() if graph.degree[ v ]  == k//2 ]  #< len(genset) ] 
    for (u,v) in parttion(graph, leaves):
        graph = nx.identified_nodes(graph, u, v , self_loops = False )  
    return graph




def intersection(plaq1, plaq2, graph):
    def edges_pla(pla):
        for j in range(len(plaq1)):
            pass# return plaq[


    #for u,v in zip( plaq1, plaq1  
            

def union_of_cayleys():

    G2 = nx.Graph()
    gf = ga.GF(107)
    genset = [ ] 
    #for i in range(2): 
    _genset = [ gf(i) for i in  [   62, 99, 49, 89   ] ]
    for g in _genset:
        print(g)
        
    DEPTH = 107
    for e,j in sampleDFS_PSL(gf, _genset , gf(1), depth = DEPTH, color = set()):
        G2.add_edge(*e)
                    
    name = "union-plt-{0}-{1}".format(1, dtime.now())
    save_graph(G2, "./img/union/{0}.{1}".format(name, "svg")) 
    pkl.dump(G2, open("./pkl/{0}.{1}".format(name, "pkl"), "bw"))  
    return G2, genset

def sample_plaquettes():
    G2, genset = union_of_cayleys()

    

    name = "closed-plt-{0}-{1}".format(i, dtime.now())
    save_graph(G_temp, "./img/closed/{0}.{1}".format(name, "svg")) 
    pkl.dump(G_temp, open("./pkl/{0}.{1}".format(name, "pkl"), "bw"))  

def sample_and_close():
    for i in range(5): 
        gf = ga.GF(401)
        p = gf( 3)

        down = np.zeros((4,4), dtype = np.int )
        down[2][2] = 1
        down[3][3] = 1

        up = np.zeros((4,4), dtype = np.int )
        up[0][0] = 1
        up[1][1] = 1

        funcs = { 
            0 : lambda g : np.kron(gf([[1,0],[0,0]] ), g) + gf(down),  
            1 : lambda g : np.kron(gf([[0,0],[0,1]] ), g) + gf(up)
        }
        flag_pass = False
        while not flag_pass:
            genset = [[funcs[j](g)  for g in randomBase(gf, gens= 4) ] for j in range(2) ]
            for Ag in genset:
                for g in Ag:
                    print(g)
        
            DEPTH = 3
            G2 = nx.Graph()
            for e,j in sampleDFS_squareComplex(gf, genset , np.kron(gf([[1,0],[0,1]]), gf([[1,0],[0,1]]), ), depth = DEPTH, color = set()):
                G2.add_edge(*e, label=j)
            for _ in range(3):
                if flag_pass:
                    break
                G_temp =  contracted_leaves_k( deepcopy(G2), genset,  k = 4)
                G_temp =  contracted_leaves_k( deepcopy(G_temp), genset,  k = 4)
                G_temp =  contracted_leaves_k( deepcopy(G_temp), genset,  k = 4)
                G_temp =  contracted_leaves_k( deepcopy(G_temp), genset,  k = 4)
                if True or len([ v for v in G_temp.nodes() if G_temp.degree[v]  < 3   ]) <= 1 : 
                    flag_pass  = True
                    name = "closed-plt-{0}-{1}".format(i, dtime.now())
                    save_graph(G_temp, "./img/closed/{0}.{1}".format(name, "svg")) 
                    pkl.dump(G_temp, open("./pkl/{0}.{1}".format(name, "pkl"), "bw"))  
                    print("@")
if __name__ == "__main__" :
   
    union_of_cayleys()
    exit(0)
    sample_and_close()

    gf = ga.GF(401)
    p = gf( 5)

    Groundtruth = [ gf([[2, 5], [12, 11]]), gf([ [2, 1],[8, 11]]),gf([ [2, 12],[5, 11]]),gf([ [2, 8],[1, 11]])]

    down = np.zeros((4,4), dtype = np.int )
    down[2][2] = 1
    down[3][3] = 1

    up = np.zeros((4,4), dtype = np.int )
    up[0][0] = 1
    up[1][1] = 1

    funcs = { 
            0 : lambda g : np.kron(gf([[1,0],[0,0]] ), g) + gf(down),  
            1 : lambda g : np.kron(gf([[0,0],[0,1]] ), g) + gf(up)
    }

    genset = [[funcs[j](g)  for g in randomBase(gf) ] for j in range(2) ]
    for Ag in genset:
        for g in Ag:
            print(g)

    DEPTH = 2 
    G2 = nx.Graph()
    for e,j in sampleDFS_squareComplex(gf, genset , np.kron(gf([[1,0],[0,1]]), gf([[1,0],[0,1]]), ), depth = DEPTH):
        G2.add_edge(*e, label=j)
    leaves = [ v for v in G2.nodes() if G2.degree[v] < len(genset) ] 
    shuffle(leaves)
    for j in range(len(leaves)//2):
        G2 = nx.contracted_nodes(G2, leaves[j], leaves[j+ (len(leaves)+1)//2 ], self_loops = False)  
    save_graph(G2, "plt.svg")

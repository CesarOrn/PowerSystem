import networkx as nx
import numpy as np


def getNetwork(name):

    if(name=="symmetric"):
        G=nx.Graph([(0,2),(0,4),(3,1),(1,5),(3,2),(2,4),(6,2),(2,7),(8,2),(3,5),(3,6),(7,3),(3,8),(5,4),(9,4),(4,10),(4,11),(5,9),(5,10),(5,11),(6,7),(6,8),(7,8),(9,10),(9,11),(10,11)])
        P=[1,1,1,1,1,1,-1,-1,-1,-1,-1,-1]

    elif(name=="smallNetwork"):
        G=nx.Graph([(1,0),(1,1),(1,2),(2,1),(2,2)])
        P=[2,-8,6]

    elif(name=="star"):
        G=nx.Graph([(0,1),(0,2),(0,3),(1,4),(4,5),(5,1),(5,6),(6,2),(7,2),(7,6),(7,8),(3,8),(8,9),(9,3),(4,9)])
        P= [3,1,1,1,-1,-1,-1,-1,-1,-1]

    elif(name=="starCom"):

        G=nx.MultiDiGraph()
        G.nodes(data=True)
        G.add_edge(0, 1, label=str(3) )
        G.add_edge(1, 0, label=str(2) )
        G.add_edge(1, 2, label=str(1) )
        G.add_edge(2, 1, label=str(2) )
        G.add_edge(2, 2, label=str(1) )
        W=[3,1,2,1,2]
        P= [3,1,-1]

    elif(name=="bottleNetwork"):
        G=nx.Graph([(1,0),(0,2),(1,2),(1,3),(1,4),(2,3),(5,4),(6,4),(6,5)])
        P= [3,-1,-1,1,-1,-2,1]

    elif(name=="bottleNetworkBroken"):

        G=nx.Graph([(1,0),(0,2),(1,2),(1,3),(2,3),(5,4),(6,4),(6,5)])
        P= [3,-1,-1,1,-1,-2,1]

    else:
        G=nx.Graph([(1,2),(0,1)])
        P= [2,-3,1]

    #'overlap' : 'false'
    G.graph['graph']={'rankdir':'TD','pad':'1,1', 'size' :'7.75,10.25!','dpi':'300'}
    G.graph['node']={'style':'filled','fixedsize': 'true','width': '0.5','height': '0.5'}
    G.graph['edge']={'arrowsize':'1.2','penwidth':'0.8'}

    return G,P

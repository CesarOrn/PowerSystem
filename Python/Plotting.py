
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import networkx as nx
from networkx.drawing.nx_agraph import to_agraph
import numpy as np
import csv
import matplotlib.cm as cm
import colorsys
from network import *
import graphviz
import math

#!/usr/bin/env python
def clamp(num, min_value, max_value):
   return max(min(num, max_value), min_value)

def trans(pos,amount=150):
    NewPos={}
    for (node, (x,y)) in pos.items():
        print(amount)
        theta = (amount/180)*(math.pi)
        print(theta)
        Rotate = np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]])
        X=float(x)
        Y=float(y)

        point =np.zeros((2,1))
        point[0,0]=X
        point[1,0]=Y
        out= np.dot(Rotate,point)
        NewPos[node]=[out[0,0],out[1,0]]
    return NewPos

networkData=getNetwork("starCom")



print('hello')
P=networkData[1]
#P=[4,-1,-2,2,-2,-1]
max_value = np.max(P)
GR=networkData[0]
#GR= nx.Graph([(0,1),(1,2),(0,2),(0,4),(3,2),(3,4),(5,4)])
#G = nx.adjacency_matrix(GR)

pos=nx.kamada_kawai_layout(GR)

pos=trans(pos,150)


scale=2.1
for (node, (x,y)) in pos.items():
    if (node ==1 or node ==2 or node ==3):
        #GR.nodes[node]['pos']=str((scale/2)*x)+","+str((scale/2)*y)+"!"
        GR.nodes[node]['pos']=str(scale*x)+","+str(scale*y)+"!"
    else:
        GR.nodes[node]['pos']=str(scale*x)+","+str(scale*y)+"!"


GR=networkData[0]
print(networkData[0].nodes[0])


#plot vertex consumer
loopSize=nx.number_of_nodes(GR)
print('size:',loopSize)
nodeSize=1500
fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(6,5))
ax = axes


for i in  range(loopSize):

    if(P[i]>=0):
        #Producer
        print('lookin at %d consumer',i)
        scale=abs(P[i]/max_value)
        print('scale is',scale)
        if i==0:
            c='g'
            #GR.nodes[i]['fillcolor']="#00bfbf"
            GR.nodes[i]['fillcolor']="#00ff00"
        else:
            GR.nodes[i]['fillcolor']="#00ff00"
            c='g'
        #node_color=colorsys.hls_to_rgb(0.30, clamp(scale, 0.0, 0.7), 1)
        nx.draw_networkx_nodes(GR,pos,nodelist=[i], node_color=c,node_size=nodeSize, alpha=1.0,node_shape='s')
        GR.nodes[i]['shape']="square"
    else:
        #Consumer
        print('lookin at %d producer',i)
        scale=abs(P[i]/max_value)
        #node_color=colorsys.hls_to_rgb(0.0, clamp(scale, 0.0, 0.9), 1)
        nx.draw_networkx_nodes(GR,pos,nodelist=[i],node_color='r', node_size=nodeSize, alpha=1.0,node_shape='o')
        GR.nodes[i]['fillcolor']="#ff0000"
        GR.nodes[i]['shape']="circle"




#plot edges
nx.draw_networkx_edges(GR,pos,width=1.25,alpha=1.0, arrows=True)


#Label
labels = {}
for i in range(loopSize):
    j=i+1
    labels[i]=(str(j))
    GR.nodes[i]['label']=str(j);


label={ i:j for i,j in labels.items() if i in pos}

nx.draw_networkx_labels(GR, pos, label, font_size=14)
#ax[0].text(0.5, -0.1, '(a)', transform=ax[0].transAxes,fontsize=16, fontweight='regular', va='bottom', ha='center')
ax.axis('off')

A = nx.drawing.nx_agraph.to_agraph(GR)
print(GR)
print(A)
A.layout()
print(A)
A.draw('multi.png')



X = []
Y = [[],[],[],[],[],[],[],[],[],[]]
TR = [[],[],[],[],[],[]]

with open('InitQPert.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    next(plots)
    for row in plots:
        X.append(float(row[0]))
        Y[0].append(float(row[1]))
        Y[1].append(float(row[2]))
        Y[2].append(float(row[3]))
        #Y[3].append(float(row[4]))
        #Y[4].append(float(row[5]))
        #Y[5].append(float(row[6]))
        #Y[6].append(float(row[7]))
        #Y[7].append(float(row[8]))
        #Y[8].append(float(row[9]))
        #Y[9].append(float(row[10]))

with open('trig.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    next(plots)
    for row in plots:
        TR[0].append(float(row[1]))
        TR[1].append(float(row[2]))
        TR[2].append(float(row[3]))
        TR[3].append(float(row[4]))
        TR[4].append(float(row[5]))
        TR[5].append(float(row[6]))
#color
#cmap = plt.get_cmap('tab10')
cmap = plt.get_cmap('gist_rainbow')
colors = [cmap(i) for i in np.linspace(0, 1.0, 3)]

plt.figure(2,figsize=(7,5))


for i, color in enumerate(colors,start=0):
    #plt.plot(X, Y[i],color=color)
    plt.plot(X, Y[i])

#for i, color in enumerate(colors,start=0):
    #plt.plot(X,TR[i], linestyle=':', color='gray')


#plt.grid(axis='x', color='0.85')
plt.grid(axis='y', color='0.75')
plt.ylim([-1,3])
ax = plt.gca()
ax.set_facecolor((0.9,0.9,0.9))
plt.xlabel('Time (S)', fontsize=14)
plt.ylabel('$\eta$', fontsize=14)
#plt.ylabel('$\theta$', fontsize=14)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
#plt.show()


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import csv
import matplotlib.cm as cm
import colorsys
from network import *

def clamp(num, min_value, max_value):
   return max(min(num, max_value), min_value)


networkData=getNetwork("smallNetwork")

print('hello')
P=networkData[1]
#P=[4,-1,-2,2,-2,-1]
max_value = np.max(P)
GR=networkData[0]
#GR= nx.Graph([(0,1),(1,2),(0,2),(0,4),(3,2),(3,4),(5,4)])
#G = nx.adjacency_matrix(GR)



pos=nx.kamada_kawai_layout(GR)

#plot vertex consumer
loopSize=nx.number_of_nodes(GR)
print('size:',loopSize)
nodeSize=1000
plt.figure(1)
for i in  range(loopSize):

    if(P[i]>=0):
        #Producer
        print('lookin at %d consumer',i)
        scale=abs(P[i]/max_value)
        print('scale is',scale)
        #node_color=colorsys.hls_to_rgb(0.30, clamp(scale, 0.0, 0.7), 1)
        nx.draw_networkx_nodes(GR,pos,nodelist=[i], node_color='g',node_size=nodeSize, alpha=1.0,node_shape='s')
    else:
        #Consumer
        print('lookin at %d producer',i)
        scale=abs(P[i]/max_value)
        #node_color=colorsys.hls_to_rgb(0.0, clamp(scale, 0.0, 0.9), 1)
        nx.draw_networkx_nodes(GR,pos,nodelist=[i],node_color='r', node_size=nodeSize, alpha=1.0,node_shape='o')




#plot edges
nx.draw_networkx_edges(GR,pos,width=1.25,alpha=1.0)


#Label
labels = {}
for i in range(loopSize):
    j=i+1
    labels[i]=(r"$\theta_{"+str(j)+"}$")


label={ i:j for i,j in labels.items() if i in pos}

nx.draw_networkx_labels(GR, pos, label, font_size=10)

#show
plt.axis("off")




X = []
Y = [[],[],[],[],[],[]]
TR = [[],[],[],[],[],[]]

with open('plot.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    next(plots)
    for row in plots:
        X.append(float(row[0]))
        Y[0].append(float(row[1]))
        Y[1].append(float(row[2]))
        Y[2].append(float(row[3]))
        Y[3].append(float(row[4]))
        Y[4].append(float(row[5]))
        Y[5].append(float(row[6]))

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
cmap = plt.get_cmap('gist_rainbow')
colors = [cmap(i) for i in np.linspace(0, 0.5, 6)]


plt.figure(2)
for i, color in enumerate(colors,start=0):
    plt.plot(X, Y[i])

#for i, color in enumerate(colors,start=0):
    #plt.plot(X,TR[i], linestyle=':', color='gray')


plt.grid(axis='x', color='0.85')
#plt.grid(axis='y', color='0.85')
plt.ylim([-3,3])
ax = plt.gca()
ax.set_facecolor((0.9,0.9,0.9))
plt.xlabel('Time(S)', fontsize=10)
plt.ylabel('Etas', fontsize=10)
plt.show()

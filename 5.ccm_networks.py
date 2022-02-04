#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 10:39:34 2021

@author: valentin
"""

from numpy import array
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import os

labelist = ["D","E","O","T","C","S","Sf"] #list of variables

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

fig = plt.figure(figsize=(10,13)) # x inch * y inch

coord = [(0, 0),(0, 1),(0, 2),(1, 0),(1, 1),(1, 2), #plot coordinates
         (2, 0),(2, 1),(2, 2),(3, 0),(3, 1),(3, 2)]

for i in range(9): #for each taxonomic dataset (folder)

    #load ccm tables
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),
             "datasets","taxonomic_databases",namefile[i]))

    with open("rho_matrix.csv", "r") as tablefile:
        tablerho = [[x for x in ln.split()] for ln in tablefile]
        
    with open("pval_matrix.csv", "r") as tablefile:
        tablepval = [[x for x in ln.split()] for ln in tablefile]
        
    #load subplot
    ax4 = fig.add_subplot(4, 3, i+1)
    ax4.set_title(namefile[i])

    G = nx.MultiDiGraph()
    G.add_nodes_from(labelist)
    pos = nx.circular_layout(G)
    
    #add arrows and nodes
    edge_colors = []
    for row in range(1,8):
        for col in range(1,4):
            
            #display arrow if p-value<0.05
            if tablerho[row][col] != "NA" and float(tablepval[row][col]) <= 0.05: 
                G.add_edge(labelist[row-1], labelist[col-1], 
                           weight=float(tablerho[row][col]))
                edge_colors.append(float(tablerho[row][col]))

    nodes = nx.draw_networkx_nodes(G, pos, node_size=400, 
            node_color=["#a9a9a9","#ff335e", "#3399ff", "#33ff86","#efeded", 
                        "#d9d9d9", "#a9a9a9"], linewidths=1,edgecolors="black")
    
    pos = {'Sf': array([1.00000000e+00, 8.51494896e-09]),
           'E': array([0.62348981, 0.78183148]),
           'D': array([-0.22252091,  0.97492788]),
           'O': array([-0.90096878,  0.43388381]),
           'T': array([-0.90096878, -0.43388373]),
           'C': array([-0.22252097, -0.97492786]),
           'S': array([ 0.62348963, -0.78183159])}
    
    nx.draw_networkx_labels(G, pos, with_labels = True,font_weight='bold')
    
    maxcolor = 0.7
    
    edges = nx.draw_networkx_edges(
        G,
        pos,
        node_size=400,
        arrowstyle="-|>",
        arrowsize=10,
        edge_color=edge_colors,
        edge_cmap=plt.cm.Blues,
        edge_vmin=0.0,
        edge_vmax=maxcolor,
        width=2,
        connectionstyle='arc3, rad = 0.1')
    
    try:
        pc = mpl.collections.PatchCollection(edges, cmap=plt.cm.Blues)
    except:
        pass

    pc.set_array(edge_colors)
    ax = plt.gca()
    ax.set_axis_off()

pc.set_clim(vmin=0, vmax=maxcolor)
plt.colorbar(pc)    
plt.savefig("ccm.pdf")
plt.show()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### V. Rineau, J, Smyčka, D. Storch, Diversity-dependence is an ubiquitous phenomenon across phanerozoic oceans.
### DATA 1.
###
### Data 1: script 5 - 5.ccm_networks.py
###
### Reproduces the Figure 1 of the paper following CCM analyses computed in script 4.

# Load packages
from numpy import array
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import os

labelist = ["D","E","O","T","C","S","Sf"] #list of variables

namefile = ["bivalvia","brachipoda","scleractinia", #list of folders
            "gasteropoda","metazoa","prymnesiophycae",
            "foraminifera","coccolithophoridae","radiolaria"]

# Load figure
fig = plt.figure(figsize=(10,13)) # x inch * y inch

coord = [(0, 0),(0, 1),(0, 2),(1, 0),(1, 1),(1, 2), 
         (2, 0),(2, 1),(2, 2),(3, 0),(3, 1),(3, 2)]

main_folder = os.path.dirname(os.path.abspath(__file__))

# For each taxonomic dataset
for i in range(9): 

    # Load convergent-cross maping results
    os.chdir(os.path.join(main_folder,
             "datasets","taxonomic_databases",namefile[i]))

    with open("rho_matrix.csv", "r") as tablefile:
        tablerho = [[x for x in ln.split()] for ln in tablefile]
        
    with open("pval_matrix.csv", "r") as tablefile:
        tablepval = [[x for x in ln.split()] for ln in tablefile]
        
    # Load subplot for specific taxaset
    ax4 = fig.add_subplot(4, 3, i+1)
    ax4.set_title(namefile[i])

    G = nx.MultiDiGraph()
    G.add_nodes_from(labelist)
    pos = nx.circular_layout(G)
    
    # Add arrows and nodes
    edge_colors = []
    for row in range(1,8):
        for col in range(1,4):
            
            # Display arrow if p-value < 0.05
            if tablerho[row][col] != "NA" and float(tablepval[row][col]) <= 0.05: 
                G.add_edge(labelist[row-1], labelist[col-1], 
                           weight=float(tablerho[row][col]))
                edge_colors.append(float(tablerho[row][col]))

    # Positionning elements / graphical settings
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
    
    nx.draw_networkx_labels(G, pos, font_weight='bold')
    
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

# Save and show
os.chdir(os.path.join(main_folder,"datasets","taxonomic_databases"))
plt.savefig("ccm.pdf")
plt.show()


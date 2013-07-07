'''
Created on May 27, 2012

@author: eskin3
'''

import LDModel
import networkx as nx
import matplotlib.pyplot as plt
import string
from networkx import graphviz_layout

def draw_HMM(h, anc = 0, start_level = 0, level_num = None):
    G=nx.DiGraph()
    G.rtt = {}
    (edges, edge_weights) = h.get_node_edges(anc)
    nodes = h.get_layer_node_nums(anc)
    ems_probs = h.get_node_ems_probs(anc)
    snp_num = h.get_snp_num()
    if level_num == None:
        level_num = snp_num
    if start_level + level_num > snp_num:
        level_num = snp_num - start_level
    #edges = edges[start_level:start_level+level_num]
    #nodes = nodes[start_level:start_level+level_num]
    for level_idx in range(start_level, start_level + level_num):
        print "adding level: " + str(level_idx) + " to the graph"
        for node_idx in range(nodes[level_idx]):
            print "  adding node: " + str(node_idx) + " to the graph"
            new_node = str(level_idx) + "." + str(node_idx)
            G.add_node(new_node)
            if ems_probs[level_idx][node_idx][0] > 0.5:
                G.rtt[new_node] = 0.3
            else:
                G.rtt[new_node] = 0.6
            if level_idx < start_level + level_num - 1:
                for other_node_idx in range(len(edges[level_idx][node_idx])):
                    other_node = edges[level_idx][node_idx][other_node_idx]
                    edge = (str(level_idx) + "." + str(node_idx),str(level_idx+1) + "." + str(other_node))
                    G.add_edge(str(level_idx) + "." + str(node_idx),str(level_idx+1) + "." + str(other_node), {'prob':edge_weights[level_idx][node_idx][other_node_idx]})
                    other_node_idx += 1
    
    pos=nx.graphviz_layout(G,prog="dot",root="0.0")
    min_x = min([x[0] for x in pos.values()])
    min_y = min([x[1] for x in pos.values()])
    max_x = max([x[0] for x in pos.values()])
    max_y = max([x[0] for x in pos.values()])
    for key in pos.keys():
        level = int(string.split(key,".")[0])
        node_num = int(string.split(key,".")[1])
        print "level: " + str(level) 
        pos[key] = ((max_x - min_x) / len(nodes) * level, (max_y - min_y) / max(nodes) * node_num)
    nx.draw(G, pos, node_color=[G.rtt[v] for v in G], alpha=0.5) #node_color=[G.rtt[v] for v in G], node_size=15
    edge_labels=dict([((u,v,),"{0:.2f}".format(d['prob'])) for u,v,d in G.edges(data=True)])
    nx.draw_networkx_edge_labels(G,pos,alpha=0.5,edge_labels=edge_labels) 
    plt.show()
    print "a"
        
def draw_old_graph():
    pass
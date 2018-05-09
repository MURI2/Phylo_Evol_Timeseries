from __future__ import division
import os, math, numbers, itertools, re
from itertools import permutations
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")


def draw_graph(graph, labels=None, graph_layout='spring',
               node_size=50, node_color='blue', node_alpha=0.3,
               node_text_size=3,
               edge_color='blue', edge_alpha=0.3, edge_tickness=0.4,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # create networkx graph
    G = nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size,
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels,
                                 label_pos=edge_text_pos, font_size = 1)

    # show graph
    #plt.show()
    plt.savefig(mydir + 'figs/metabolic_networks/test_network.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plot_network():
    IN_path = mydir + 'data/reference_assemblies_task2/MAPLE/MAPLE_modules/B_MAPLE_modules/B_KO_M_protein.txt'
    IN = pd.read_csv(IN_path, sep = '\t')
    to_iter = IN[(IN.Pathway_ID.notnull()) & (IN.protein_id.notnull())]
    to_iter_group = to_iter.groupby(['Pathway_ID'])
    #test =  to_iter_group.get_group('M00793')
    graph = []
    for name, group in to_iter_group:
        kegg = group.KEGG_Orthology.values
        len_legg = len(kegg)
        for i in range(1, len_legg):
            for j in range(0, i):
                graph.append((kegg[i], kegg[j]))

    node_dict = {}
    #labels = map(int, range(65, 65+len(graph)))
    #draw_graph(graph)

def merge_pathways():
    strains = ['B', 'C', 'D', 'F', 'J', 'P']
    # go through, make dict of dict[pathway] = list of genes
    # then make one large pathway 
    print "work on this"




#plot_network()
plot_network()

from datetime import datetime
import numpy as np
import igraph as ig
from nexus import NexusReader

PLOT_NUM = 1
def myplot(obj, *args, **kwargs):
    global PLOT_NUM
    ig.plot(obj, target="cairo_dumps/"+datetime.now().strftime('%d-%H-%M-%S')+"__{0}.png".format(PLOT_NUM),
                    orientation="right-left", *args, **kwargs)
    PLOT_NUM += 1

# returns dictionary where values are lists of strings
def nex2char_dict(nex_file):
    n = NexusReader()
    n.read_file(nex_file)
    return n.data.matrix # a dict

def char_dict2dist_mat(char_dict):
    def dist(a,b):
        pass

    n = len(char_dict)
    distances = np.zeros((n,n))
    keys = char_dict.keys()
    for ix1, key1 in enumerate(keys):
        for ix2, key2 in enumerate(keys):
            distances[ix1, ix2] = dist(char_dict[key1], char_dict[key2])

    return distances

# ig.Graph.Adjacency(DistMatrix, mode = 'DIRECTED') won't work with weights
def build_lang_graph(dist_matrix, labels):
    n = len(labels)
    dist_ray = np.array(dist_matrix)
    assert n == dist_ray.shape[0] == dist_ray.shape[1]
    assert np.allclose(dist_ray, dist_ray.T, atol=0.01)
    L = ig.Graph.Full(n)
    L.vs['label'] = labels
    for i in range(n):
        for j in range(i+1, n):
            L.es[L.get_eid(i,j)]["weight"] = dist_ray[i,j]
    return L

def plot_comms(clustering, *args, **kwargs):
    myplot(clustering, vertex_color=clustering.membership, *args, **kwargs)

"""
g = Graph()
g.add_vertices(3) #adds 3 edges
g.add_edges([(0,1), (1,2)])
g.vs["label"] = ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"] #labels
g.es[g.get_eid(2,3)]["is_formal"] = True #adds a label two the specific edge 2~3

"""

"""
#methods that return a dengrogram -- see https://igraph.org/python/doc/igraph.clustering.VertexDendrogram-class.html
community_fastgreedy(self, weights=None)
community_leading_eigenvector_naive(clusters=None, return_merges=False)
community_edge_betweenness(self, clusters=None, directed=True, weights=None)
community_walktrap(self, weights=None, steps=4)
FURTHER, we can adapt: community_multilevel(self, weights=None, return_levels=False)
"""

karate = ig.load("karate.gml")
layout = karate.layout("kk")
cl = karate.community_leading_eigenvector(clusters=2)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector(clusters=3)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector(clusters=4)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector()
plot_comms(cl, layout=layout)

#cl = karate.community_leading_eigenvector_naive()


from turk_given_dist import Turk_dist_mat, Turk_labels
g = build_lang_graph(Turk_dist_mat, Turk_labels)
layout = g.layout("kk")
cl = g.community_leading_eigenvector(weights=g.es['weight'])
plot_comms(cl, layout=layout)

"""
myplot(g, layout = layout)

cl = g.community_fastgreedy(weights=g.es['weight'])
myplot(cl, layout=layout)
plot_comms(cl.as_clustering())
#myplot(cl, layout=layout, vertex_label=cl.as_clustering().membership)"""

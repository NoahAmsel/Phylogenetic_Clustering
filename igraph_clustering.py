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

def plot_comms(clustering, *args, **kwargs):
    myplot(clustering, vertex_color=clustering.membership, *args, **kwargs)

# returns dictionary where values are lists of strings
def nex2char_dict(nex_file, missing="-"):
    n = NexusReader()
    n.read_file(nex_file)
    char_dict = {}
    for label, chars in n.data.matrix.iteritems(): #a dict
        char_dict[label.encode('ascii','replace')] = np.array([np.nan if x==missing else float(x) for x in chars])
    return char_dict

def hamming_sim(a,b):
    #assumes we're getting np.arrays
    assert a.shape == b.shape
    shared_features = (a==b).sum()
    present_features = (~(np.isnan(a)|np.isnan(b))).sum()
    return float(shared_features)/present_features

def inverse_hamming_dist(a,b):
    #assumes we're getting np.arrays
    assert a.shape == b.shape
    present_features = (~(np.isnan(a)|np.isnan(b))).sum()
    unshared_features = (a!=b).sum() - (len(a) - present_features) #np.nan is unequal to everything
    if unshared_features==0:
        return float(present_features)
    else:
        return float(present_features)/unshared_features

def char_dict2sim_mat(char_dict, sim_fun):
    n = len(char_dict)
    similarities = np.zeros((n,n))
    keys = char_dict.keys()
    for ix1, key1 in enumerate(keys):
        for ix2, key2 in enumerate(keys):
            similarities[ix1, ix2] = sim_fun(char_dict[key1], char_dict[key2])
    return similarities

# ig.Graph.Adjacency(DistMatrix, mode = 'DIRECTED') won't work with weights
def build_lang_graph(sim_matrix, labels): #sim_matrix is SIMILARITY scores not distances
    n = len(labels)
    sim_ray = np.array(sim_matrix)
    assert n == sim_ray.shape[0] == sim_ray.shape[1]
    assert np.allclose(sim_ray, sim_ray.T, atol=0.001)
    L = ig.Graph.Full(n)
    L.vs['label'] = labels
    for i in range(n):
        for j in range(i+1, n):
            L.es[L.get_eid(i,j)]['weight'] = sim_ray[i,j]
    L.es['width'] = L.es['weight']
    return L

def nexus2lang_graph(nex_file, sim_fun, missing="-"):
    chars = nex2char_dict(nex_file, missing=missing)
    sims = char_dict2sim_mat(chars, sim_fun)
    return build_lang_graph(sims, chars.keys())

def dendro2nexus(dendro, file):
    with open(file, 'w') as f:
        graph = dendro.as_clustering().graph
        f.write('#NEXUS\n')
        f.write('Begin trees;\n')
        if 'label' in graph.vs.attribute_names():
            f.write('\ttranslate\n')
            for ix, v in enumerate(graph.vs):
                if ix==len(graph.vs)-1:
                    # semicolon
                    f.write('\t\t{0}\t{1};\n'.format(ix, v['label']))
                else:
                    f.write('\t\t{0}\t{1},\n'.format(ix, v['label']))
        f.write('\t\ttree Dendrogram = {0}\n'.format(dendro.format()))
        f.write('end;')

def nested2nexus(nested, labels, file):
    def iterative_newkirk(nested):
        if not hasattr(nested, '__iter__'):
            return str(int(nested))
        else:
            return "(" + ",".join([iterative_newkirk(x) for x in nested]) + ")"

    with open(file, 'w') as f:
        f.write('#NEXUS\n')
        f.write('Begin trees;\n')
        f.write('\ttranslate\n')
        for ix, v in enumerate(labels):
            if ix==len(labels)-1:
                # semicolon
                f.write('\t\t{0}\t{1};\n'.format(ix, v))
            else:
                f.write('\t\t{0}\t{1},\n'.format(ix, v))
        f.write('\t\ttree Dendrogram = {0};\n'.format(iterative_newkirk(nested)))
        f.write('end;')

def nested2dendro(nested, graph):
    nested2dendro.merges = []
    nested2dendro.next_vert = len(graph.vs)

    def helper(nested):
        if hasattr(nested, '__iter__'):
            assert len(nested) == 2
            nested2dendro.merges.append((helper(nested[0]),helper(nested[1])))
            nested2dendro.next_vert += 1
            return nested2dendro.next_vert-1
        else:
            return nested

    helper(nested)
    return ig.VertexDendrogram(graph, nested2dendro.merges)

def iterative_community_binary(graph):
    if len(graph.vs)==1:
        return graph.vs['id'][0]

    if 'id' not in graph.vs.attributes():
        graph.vs['id'] = list(range(1,1+len(graph.vs)))
    w = 'weight' if 'weight' in graph.es.attributes() else None
    cl = graph.community_leading_eigenvector(weights=w, clusters=2)
    if len(cl.subgraphs())==1:
        cl = graph.community_fastgreedy(weights=w).as_clustering(n=2)
    return [iterative_community_binary(sub) for sub in cl.subgraphs()]

def iterative_community(graph, clusters=None):
    if 'id' not in graph.vs.attributes():
        graph.vs['id'] = list(range(1,1+len(graph.vs)))
    w = 'weight' if 'weight' in graph.es.attributes() else None
    cl = graph.community_leading_eigenvector(weights=w, clusters=clusters)
    if len(cl.subgraphs())==1:
        if len(graph.vs)==1:
            return graph.vs['id'][0]
        else:
            return graph.vs['id']
    else:
        return [iterative_community(sub, clusters=clusters) for sub in cl.subgraphs()]





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

"""karate = ig.load("karate.gml")
layout = karate.layout("kk")
cl = karate.community_leading_eigenvector(clusters=2)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector(clusters=3)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector(clusters=4)
plot_comms(cl, layout=layout)
cl = karate.community_leading_eigenvector()
plot_comms(cl, layout=layout)"""

#cl = karate.community_leading_eigenvector_naive()


"""from turk_given_dist import Turk_dist_mat, Turk_labels
g = build_lang_graph(Turk_dist_mat, Turk_labels)
layout = g.layout("kk")
cl = g.community_leading_eigenvector(weights=g.es['weight'])
plot_comms(cl, layout=layout)"""

"""
myplot(g, layout = layout)

cl = g.community_fastgreedy(weights=g.es['weight'])
myplot(cl, layout=layout)
plot_comms(cl.as_clustering())
#myplot(cl, layout=layout, vertex_label=cl.as_clustering().membership)
"""

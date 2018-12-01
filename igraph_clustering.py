from datetime import datetime
from collections import defaultdict
import numpy as np
import igraph as ig
from nexus import NexusReader

def myplot(obj, *args, **kwargs):
    ig.plot(obj, target="cairo_dumps/"+datetime.now().strftime('%d-%H-%M-%S')+"__{0}.png".format(myplot.PLOT_NUM),
                    orientation="right-left", *args, **kwargs)
    myplot.PLOT_NUM += 1
myplot.PLOT_NUM = 1

def plot_comms(clustering, *args, **kwargs):
    myplot(clustering, vertex_color=clustering.membership, *args, **kwargs)

karate = ig.load("karate.gml")

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
    """
    DEPRECATED
    """
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

def labelednested2nexus(nested, labels, file, ix_shift=0):
    """
    DEPRECATED
    """
    def iterative_newkirk(nested):
        if not hasattr(nested, '__iter__'):
            return str(int(nested+ix_shift))
        else:
            return "(" + ",".join([iterative_newkirk(x) for x in nested]) + ")"

    with open(file, 'w') as f:
        f.write('#NEXUS\n')
        f.write('Begin trees;\n')
        f.write('\ttranslate\n')
        for ix, v in enumerate(labels):
            if ix==len(labels)-1:
                # semicolon
                f.write('\t\t{0}\t{1};\n'.format(ix+ix_shift, v))
            else:
                f.write('\t\t{0}\t{1},\n'.format(ix+ix_shift, v))
        f.write('\t\ttree Dendrogram = {0};\n'.format(iterative_newkirk(nested)))
        f.write('end;')



def nested2dendro(nested, graph):
    nested2dendro.merges = []
    nested2dendro.next_vert = graph.vcount()

    def helper(nested):
        if hasattr(nested, '__iter__'):
            assert len(nested) == 2
            nested2dendro.merges.append((helper(nested[0]),helper(nested[1])))
            nested2dendro.next_vert += 1
            return nested2dendro.next_vert-1
        else:
            return int(nested)

    helper(nested)
    return ig.VertexDendrogram(graph, nested2dendro.merges)

def fast_greedy_wrapper(graph, weights=None):
    dendro = graph.community_fastgreedy(weights=weights)
    return list(range(graph.vcount())), dendro.merges

def edge_betweenness_wrapper(graph, weights=None):
    dendro = graph.community_edge_betweenness(directed=False, weights=weights)
    return list(range(graph.vcount())), dendro.merges

def safe_edge_betweenness_wrapper(graph, weights=None):
    if graph.vcount() > 40:
        return fast_greedy_wrapper(graph, weights=weights)
    else:
        return edge_betweenness_wrapper(graph, weights=weights)

def newman_wrapper(graph, weights=None, arpack_options=None):
    """
    Same as community_leading_eigenvector method of graph,
    but modified to return the merge history
    return is tuple of membership list, merge history (where numbers are the cluster numbers)
    """

    kwds = dict(weights=weights)
    if arpack_options is not None:
        kwds["arpack_options"] = arpack_options

    cluster_list, merges, _ = ig.GraphBase.community_leading_eigenvector(graph, -1, **kwds)
    return cluster_list, merges

def newman_tree(graph, method, weights=None, labels=None, backup=None):
    """
    method and backup are functions that take a graph and an optional
    weights parameter and return a nested list
    """
    if isinstance(weights, basestring):
        weights = graph.es[weights]
    if labels == None:
        labels = [str(i) for i in range(graph.vcount())]
    if backup == None:
        backup = lambda graph, weights, labels: labels
    assert len(labels) == graph.vcount()

    # base case
    if graph.vcount() == 1:
        return labels[0]

    cluster_list, merges = method(graph, weights=weights)

    if len(merges)==0:
        # method didn't find anything to split
        # try the backup method
        return newman_tree(graph, backup, weights=weights, labels=labels, backup=None)
    else:
        # collect the clusters
        cluster2subtree = defaultdict(list)
        for vertex_ix, cluster_ix in enumerate(cluster_list):
            cluster2subtree[cluster_ix].append(vertex_ix)

        # resolve each one to a nested list
        for cluster_ix, vertex_list in cluster2subtree.iteritems():
            subgraph = graph.subgraph(vertex_list)
            cluster2subtree[cluster_ix] = newman_tree(subgraph,
                                                      method,
                                                      weights=weights,
                                                      labels=[labels[v] for v in sorted(vertex_list)],
                                                      backup=backup)

        next_inner_node_ix = max(cluster2subtree.keys())+1
        for left, right in merges:
            cluster2subtree[next_inner_node_ix] = [cluster2subtree[left], cluster2subtree[right]]
            next_inner_node_ix += 1

        return cluster2subtree[max(cluster2subtree.keys())]

"""#print newman_tree(karate, backup=lambda)
plot_comms(karate.community_leading_eigenvector(), vertex_label=range(karate.vcount()))
t = newman_tree(karate, newman_wrapper, backup=fast_greedy_wrapper)
print(t)
myplot(nested2dendro(t, karate), vertex_label=range(karate.vcount()))"""

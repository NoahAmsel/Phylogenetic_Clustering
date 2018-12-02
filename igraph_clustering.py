from datetime import datetime
from collections import defaultdict
import numpy as np
import igraph as ig
from nexus import NexusReader
from pandas import read_csv
from sklearn.cluster.bicluster import SpectralCoclustering

def taxa_list(taxa_csv):
    return list(read_csv(taxa_csv)['taxon'])

def myplot(obj, *args, **kwargs):
    ig.plot(obj, target="cairo_dumps/"+datetime.now().strftime('%d-%H-%M-%S')+"__{0}.png".format(myplot.PLOT_NUM),
                    orientation="right-left", *args, **kwargs)
    myplot.PLOT_NUM += 1
myplot.PLOT_NUM = 1

def plot_comms(clustering, *args, **kwargs):
    weight = 'weight' if 'weight' in clustering.graph.es.attributes() else None
    if "layout" not in kwargs:
        kwargs['layout'] = clustering.graph.layout_fruchterman_reingold(weights=weight)
    myplot(clustering, vertex_color=clustering.membership, *args, **kwargs)

karate = ig.load("karate.gml")

# returns dictionary where values are lists of strings
def nex2char_dict(nex_file, missing="-", taxa=None):
    n = NexusReader()
    n.read_file(nex_file)
    char_dict = {}
    for label, chars in n.data.matrix.iteritems(): #a dict
        if taxa is None or label in taxa:
            char_dict[label.encode('ascii','replace')] = np.array([np.nan if x==missing else float(x) for x in chars])
    print("Converting labels to ascii -- is this necessary?")
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
            if sim_ray[i,j] > 0:
                L.es[L.get_eid(i,j)]['weight'] = sim_ray[i,j]
            else:
                L.delete_edges(L.get_eid(i,j))
    L.es['width'] = L.es['weight']
    return L

def nexus2lang_graph(nex_file, sim_fun, missing="-", taxa=None):
    chars = nex2char_dict(nex_file, missing=missing, taxa=taxa)
    sims = char_dict2sim_mat(chars, sim_fun)
    return build_lang_graph(sims, chars.keys())

def nex2bipart_graph(nex_file, missing="-", taxa=None):
    lang2chars = nex2char_dict(nex_file, missing=missing, taxa=taxa)

    """
    n_langs = len(lang2chars)
    n_chars = len(lang2chars.values()[0])
    mat = np.array([v for _, v in lang2chars.iteritems()])
    adj = np.block([[np.zeros((n_langs, n_langs)), mat], [mat.T, np.zeros((n_chars, n_chars))]])
    bi = ig.Graph.Adjacency(adj.tolist(), ig.ADJ_MAX)
    """

    bi = ig.Graph()
    n_langs = len(lang2chars)
    n_chars = len(lang2chars.values()[0])
    lang_labels = lang2chars.keys()
    bi.add_vertices(lang_labels)
    bi.add_vertices(["char{0}".format(i) for i in range(n_chars)])
    bi.vs['label'] = bi.vs['name']
    bi.vs['type'] = [True]*n_langs + [False]*n_chars
    for lang in lang2chars.keys():
        bi.add_edges([(lang, "char{0}".format(i)) for i in np.where(lang2chars[lang])[0]])

    singleton_word = bi.vs.select(_degree = 0)
    bi.delete_vertices(singleton_word)
    return bi, lang_labels

# DO NOT DELETE THIS
ig.Graph.n_langs = lambda self: sum(self.vs['type'])
ig.Graph.n_chars = lambda self: self.vcount()-sum(self.vs['type'])

def get_bipart_mat(bi):
    mat = np.array(bi.get_adjacency(type=ig.GET_ADJACENCY_UPPER).data)
    return mat[:bi.n_langs(), bi.n_langs():]


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

def nested2nexus(nested, file):
    def to_newkirk(nested):
        if not hasattr(nested, '__iter__'):
            return nested
        else:
            return "(" + ",".join([to_newkirk(x) for x in nested]) + ")"

    with open(file, 'w') as f:
        f.write('#NEXUS\n')
        f.write('Begin trees;\n')
        f.write('\t\ttree Dendrogram = {0};\n'.format(to_newkirk(nested)))
        f.write('end;')

def nested2dendro(nested, graph, label):
    if isinstance(label, basestring):
        label2ix = {label: int(ix) for ix, label in enumerate(graph.vs[label])}

    nested2dendro.merges = []
    nested2dendro.next_vert = graph.vcount()

    def helper(nested):
        if hasattr(nested, '__iter__'):
            assert len(nested) == 2
            nested2dendro.merges.append((helper(nested[0]),helper(nested[1])))
            nested2dendro.next_vert += 1
            return nested2dendro.next_vert-1
        else:
            return label2ix[nested]

    helper(nested)
    return ig.VertexDendrogram(graph, nested2dendro.merges)

def prune_nested(nested, good_labels):
    if not hasattr(nested, '__iter__'):
        return nested if nested in good_labels else None
    else:
        pruned = []
        for x in nested:
            xp = prune_nested(x, good_labels)
            if xp is not None:
                pruned.append(xp)
        if len(pruned) == 0:
            return None
        elif len(pruned) == 1:
            return pruned[0]
        else:
            return tuple(pruned)

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

def coclustering_wrapper(graph, weights=None):
    X = get_bipart_mat(graph)
    if X.shape[0] <= 1 or X.shape[1] <= 1:
        return [], []
    elif X.shape[0] == 2:
        return [0,1], [(0,1)]
    else:
        try:
            clustering = SpectralCoclustering(n_clusters=2, random_state=0).fit(X)
        except:
            print(X)
        row_labels = list(clustering.row_labels_)
        if max(row_labels)-min(row_labels) == 0:
            return [], []
        return row_labels + list(clustering.column_labels_), [(0,1)]

def newman_wrapper(graph, weights=None, arpack_options=None):
    """
    Same as community_leading_eigenvector method of graph,
    but modified to return the merge history
    return is tuple of membership list, merge history (where numbers are the cluster numbers)
    """

    kwds = dict(weights=weights)
    if arpack_options is not None:
        kwds["arpack_options"] = arpack_options
    try:
        cluster_list, merges, _ = ig.GraphBase.community_leading_eigenvector(graph, -1, **kwds)
        return cluster_list, merges
    except InternalError():
        return [], []

def iterative_clustering(graph, method, weights=None, labels=None, backup=None):
    """
    method and backup are functions that take a graph and an optional
    weights parameter and return a nested list
    """
    if weights is None or isinstance(weights, basestring):
        #weights = graph.es[weights]
        pass
    else:
        assert(False)
    if isinstance(labels, basestring):
        labels = graph.vs[labels]
    if labels is None:
        labels = [str(i) for i in range(graph.vcount())]
    assert len(labels) == graph.vcount()

    # base case
    if graph.vcount() == 1:
        return labels[0]

    # deal with orphan vertices separately bc they screw up the coclustering
    orphans = graph.vs.select(_degree_lt=1)
    if len(orphans)>0:
        cluster_list = [0]*graph.vcount()
        cluster_list[orphans[0].index] = 1
        merges = [(0,1)]
    else:
        cluster_list, merges = method(graph, weights=weights)

    if len(merges)==0:
        # method didn't find anything to split
        if backup is None:
            # leave the tree unresolved
            return labels
        else:
            # try the backup method if we have one
            return iterative_clustering(graph, backup, weights=weights, labels=labels, backup=None)
    else:
        # collect the clusters
        cluster2subtree = defaultdict(list)
        for vertex_ix, cluster_ix in enumerate(cluster_list):
            cluster2subtree[cluster_ix].append(vertex_ix)

        # resolve each one to a nested list
        for cluster_ix, vertex_list in cluster2subtree.iteritems():
            subgraph = graph.subgraph(vertex_list)
            cluster2subtree[cluster_ix] = iterative_clustering(subgraph,
                                                      method,
                                                      weights=weights,
                                                      labels=[labels[v] for v in sorted(vertex_list)],
                                                      backup=backup)

        next_inner_node_ix = max(cluster2subtree.keys())+1
        for left, right in merges:
            cluster2subtree[next_inner_node_ix] = (cluster2subtree[left], cluster2subtree[right])
            next_inner_node_ix += 1

        return cluster2subtree[max(cluster2subtree.keys())]

if __name__ == "__main__":
    t = iterative_clustering(karate, newman_wrapper, backup=fast_greedy_wrapper)
    print(t)
    myplot(nested2dendro(t, karate), vertex_label=range(karate.vcount()))
    nested2nexus(t, "cairo_dumps/huh.txt")

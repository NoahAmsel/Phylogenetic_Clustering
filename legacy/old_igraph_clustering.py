def iterative_community_binary(graph):
    if len(graph.vs)==1:
        return graph.vs['id'][0]

    if 'id' not in graph.vs.attributes():
        graph.vs['id'] = list(range(len(graph.vs)))
    w = 'weight' if 'weight' in graph.es.attributes() else None
    cl = graph.community_leading_eigenvector(weights=w, clusters=2)
    if len(cl.subgraphs())==1:

        cl = graph.community_fastgreedy(weights=w).as_clustering(n=2)
    return [iterative_community_binary(sub) for sub in cl.subgraphs()]

def iterative_community(graph, clusters=None):
    if 'id' not in graph.vs.attributes():
        graph.vs['id'] = list(range(len(graph.vs)))
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

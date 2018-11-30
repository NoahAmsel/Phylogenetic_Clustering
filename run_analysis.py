from igraph_clustering import *

turk_graph = nexus2lang_graph('data/Turkic12langs.nex', inverse_hamming_dist)
bantu = nexus2lang_graph('./data/Grollemund-et-al_Bantu-database_2015_UTF.nex', inverse_hamming_dist)

def analyze_language_graph(lang_graph, folder):
    newman_comms = lang_graph.community_leading_eigenvector(weights='weight')
    greedy_tree = lang_graph.community_fastgreedy(weights='weight')
    newman_nested_list = iterative_community_binary(lang_graph)
    newman_tree = nested2dendro(newman_nested_list,lang_graph)

    dendro2nexus(greedy_tree, folder+"/greedy_tree.nex")
    nested2nexus(newman_nested_list, lang_graph, folder+"/newman_tree.nex")

    ig.plot(newman_comms, target=folder+"/comms.png", vertex_color=newman_comms.membership, edge_width=[0.1]*len(newman_comms.graph.es))
    ig.plot(greedy_tree, target=folder+"/greedy_tree.png", orientation="right-left")
    ig.plot(newman_tree, target=folder+"/newman_tree.png", orientation="right-left")

analyze_language_graph(turk_graph, "trees/turkic")
analyze_language_graph(bantu, "trees/bantu")

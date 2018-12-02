from igraph_clustering import *

#turk_complete = nexus2lang_graph('data/Turkic12langs.nex', inverse_hamming_dist)
#bantu_complete = nexus2lang_graph('./data/Grollemund-et-al_Bantu-database_2015_UTF.nex', inverse_hamming_dist)
pama_complete = nexus2lang_graph('./data/pny/pny_coded.nex', inverse_hamming_dist, taxa=taxa_list('./data/pny/pny_taxa.nex'))

def analyze_language_graph(lang_graph, folder):


    newman_comms = lang_graph.community_leading_eigenvector(weights='weight')
    layout = lang_graph.layout_fruchterman_reingold(weights='weight')
    ig.plot(newman_comms, target=folder+"/comms.png", vertex_color=newman_comms.membership,
            layout=layout, edge_width=[0.1]*len(newman_comms.graph.es), bbox=(max(600,5*lang_graph.vcount()),max(600,5*lang_graph.vcount())))

    fastgreedy_tree = iterative_clustering(lang_graph, fast_greedy_wrapper, weights='weight', labels='label')
    nested2nexus(fastgreedy_tree, folder+"/fastgreedy.nex")

    newman_tree = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label')
    nested2nexus(newman_tree, folder+"/pure_newman.nex")

    newman_tree_betweenness = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label', backup=safe_edge_betweenness_wrapper)
    nested2nexus(newman_tree_betweenness, folder+"/newman_betweenness.nex")
    ig.plot(nested2dendro(newman_tree_betweenness, lang_graph, 'label'), bbox=(600, max(600, 10*lang_graph.vcount())),
            target=folder+"/newman_btwnness.png", orientation="right-left")

    newman_tree_fastgreedy = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label', backup=fast_greedy_wrapper)
    nested2nexus(newman_tree_fastgreedy, folder+"/newman_fastgreedy.nex")
    ig.plot(nested2dendro(newman_tree_fastgreedy, lang_graph, 'label'), bbox=(600, max(600, 10*lang_graph.vcount())),
            target=folder+"/newman_greedy.png", orientation="right-left")

#analyze_language_graph(turk_complete, "trees/turkic")
#analyze_language_graph(bantu_complete, "trees/bantu")
analyze_language_graph(pama_complete, "trees/pama_nyungan")

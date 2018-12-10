from igraph_clustering import *

def analyze_language_graph(lang_graph, folder):

    newman_comms = lang_graph.community_leading_eigenvector(weights='weight')
    layout = lang_graph.layout_fruchterman_reingold(weights='weight')
    ig.plot(newman_comms, target=folder+"/comms.png", vertex_color=newman_comms.membership,
            layout=layout, edge_width=[0.1]*len(newman_comms.graph.es), bbox=(max(600,5*lang_graph.vcount()),max(600,5*lang_graph.vcount())))

    fastgreedy_tree = iterative_clustering(lang_graph, fast_greedy_wrapper, weights='weight', labels='label')
    nested2nexus(fastgreedy_tree, folder+"/fastgreedy.nex")

    newman_tree = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label')
    nested2nexus(newman_tree, folder+"/pure_newman.nex")

    newman_nonbinary_tree = iterative_clustering(lang_graph, newman_unresolved_wrapper, weights='weight', labels='label')
    nested2nexus(newman_nonbinary_tree, folder+"/newman_nonbinary.nex")

    newman_tree_betweenness = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label', backup=safe_edge_betweenness_wrapper)
    nested2nexus(newman_tree_betweenness, folder+"/newman_betweenness.nex")
    ig.plot(nested2dendro(newman_tree_betweenness, lang_graph, 'label'), bbox=(600, max(600, 10*lang_graph.vcount())),
            target=folder+"/newman_btwnness.png", orientation="right-left")

    newman_tree_fastgreedy = iterative_clustering(lang_graph, newman_wrapper, weights='weight', labels='label', backup=fast_greedy_wrapper)
    nested2nexus(newman_tree_fastgreedy, folder+"/newman_fastgreedy.nex")
    ig.plot(nested2dendro(newman_tree_fastgreedy, lang_graph, 'label'), bbox=(600, max(600, 10*lang_graph.vcount())),
            target=folder+"/newman_greedy.png", orientation="right-left")

def analyze_lang_word_graph(bigraph, lang_labels, folder):
    only_langs = bigraph.subgraph(list(range(bigraph.n_langs())))

    newman_tree_fastgreedy = prune_nested(iterative_clustering(bigraph, newman_wrapper, labels='label', backup=fast_greedy_wrapper), lang_labels)
    nested2nexus(newman_tree_fastgreedy, folder+"/bi_newman_fastgreedy.nex")
    ig.plot(nested2dendro(newman_tree_fastgreedy, only_langs, 'label'), bbox=(600, max(600, 10*bigraph.n_langs())),
            target=folder+"/bi_newman_greedy.png", orientation="right-left")

    coclustering_tree = prune_nested(iterative_clustering(bigraph, coclustering_wrapper, labels='label', backup=None), lang_labels)
    nested2nexus(coclustering_tree, folder+"/bi_pure_coclustering.nex")

    coclustering_greedy_tree = prune_nested(iterative_clustering(bigraph, coclustering_wrapper, labels='label', backup=fast_greedy_wrapper), lang_labels)
    nested2nexus(coclustering_greedy_tree, folder+"/bi_coclustering_greedy.nex")
    ig.plot(nested2dendro(coclustering_greedy_tree, only_langs, 'label'), bbox=(600, max(600, 10*bigraph.n_langs())),
            target=folder+"/bi_coclustering_greedy.png", orientation="right-left")

if __name__ == "__main__":
    turk_complete = nexus2lang_graph('data/turkic/Turkic12langs.nex', inverse_hamming_dist)
    analyze_language_graph(turk_complete, "results/turkic")

    bantu_complete = nexus2lang_graph('./data/bantu/bantu_data.nex', inverse_hamming_dist)
    analyze_language_graph(bantu_complete, "results/bantu")

    pama_complete = nexus2lang_graph('./data/pny/pny_data.nex', inverse_hamming_dist, taxa=taxa_list('./data/pny/pny_taxa.csv'))
    analyze_language_graph(pama_complete, "results/pny")

    turk_bi, turk_lang_labels = nex2bipart_graph("data/turkic/Turkic12langs.nex")
    analyze_lang_word_graph(turk_bi, turk_lang_labels, "results/turkic/bipartite")

    # can only do this for turkic, the others are too big
    biturk_newman_comms = turk_bi.community_leading_eigenvector()
    biturk_layout = turk_bi.layout_reingold_tilford(root=list(range(turk_bi.n_langs())))
    ig.plot(biturk_newman_comms, target="results/turkic/bipartite"+"/bi_comms.png", vertex_color=biturk_newman_comms.membership,
            vertex_size=[100]*turk_bi.n_langs()+[20]*turk_bi.n_chars(),
            layout=biturk_layout, edge_width=[0.1]*len(biturk_newman_comms.graph.es), vertex_frame_width=0,
            bbox=(max(600,2*turk_bi.vcount()),max(600,turk_bi.vcount())))


    bantu_bi, bantu_lang_labels = nex2bipart_graph("data/bantu/bantu_data.nex")
    analyze_lang_word_graph(bantu_bi, bantu_lang_labels, "results/bantu/bipartite")

    pny_bi, pny_lang_labels = nex2bipart_graph("data/pny/pny_data.nex", taxa=taxa_list('./data/pny/pny_taxa.csv'))
    analyze_lang_word_graph(pny_bi, pny_lang_labels, "results/pny/bipartite")

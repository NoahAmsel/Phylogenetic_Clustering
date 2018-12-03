# Phylogenetic_Clustering



Resources to get Cairo working
https://www.cairographics.org/download/
https://stackoverflow.com/questions/12072093/python-igraph-plotting-not-available -- i followed the bottom answer
https://yoyoinwanderland.github.io/Community-Detection/

http://etetoolkit.org/docs/2.3/tutorial/tutorial_trees.html#creating-trees-from-scratch


bantu data https://github.com/D-PLACE/dplace-data/tree/master/phylogenies/grollemund_et_al2015

biclustering https://scikit-learn.org/stable/modules/biclustering.html

https://en.wikipedia.org/wiki/Quantitative_comparative_linguistics#History

PNY data: ran the perl script on the tsv she sent. also had to convert the whole file to UTF-8 encoding and
correct some mis-capitalized labels (that created mismatches with the reference tree)

numerical problem with ARPACK

Python tasks
export the trees into R
get the newman method working after
bipartite clustering
running more clustering


parameters:
-- which family or subfamily (3)
-- which graph (language or bipartite)
-- which similarity function (if lang graph)
-- which splitter:
    -- fast greedy
    -- newman (backup fast greedy)
    -- newman (backup safe edge betweenness)
    -- bipartite clustering
    -- (min cut)

from modularity_maximization import partition
from modularity_maximization.utils import get_modularity

import networkx as nx
# import igraph as ig
import matplotlib.pyplot as plt
import numpy as np


# cluster_leading_eigen

np.random.seed(10)

DistMatrix = np.random.rand(40,40)
np.fill_diagonal(DistMatrix,0)
DistMatrix = DistMatrix + DistMatrix.T # we want it symmetric


avg_value = np.mean(DistMatrix)
# DistMatrix[DistMatrix < avg_value] = 0
# DistMatrix[DistMatrix > avg_value] = 1

# # DistMatrix[0:20,20:40] = .000
# # DistMatrix[20:40,0:20] = .000

print(DistMatrix)


# THIS IS USING iGRAPH
# DistMatrix = list(DistMatrix)
# np.ndarrray.tolist(DistMatrix)
# G = ig.Graph.Adjacency(DistMatrix, mode = 'DIRECTED')


# THIS IS USING NX

G = nx.from_numpy_matrix(DistMatrix,parallel_edges=True,create_using = nx.MultiGraph())

print(nx.info(G))
# comm_dict = partition(G)
# for comm in set(comm_dict.values()):
#     print("Community %d"%comm)
#     print(str([node for node in comm_dict if comm_dict[node] == comm]))


# nx.draw(G)
# plt.show()
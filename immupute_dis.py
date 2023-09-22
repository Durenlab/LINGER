import pandas as pd
import numpy as np
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix

def immupute_dis(data, data_pca, k):
    sim = squareform(pdist(data_pca))
    sim = np.max(sim) - sim
    KK = np.sum(sim, axis=0)
    twom = np.sum(KK)
    sim_norm = sim - np.outer(KK, KK) / twom
    f = np.argsort(sim_norm, axis=0)[::-1]
    KNN_P = np.column_stack((np.repeat(np.arange(sim_norm.shape[1]), k), f[:k, :].T.flatten()))
    KNN = csr_matrix((np.ones(KNN_P.shape[0]) / k, (KNN_P[:, 1], KNN_P[:, 0])), shape=sim_norm.shape)
    X = KNN.T.dot(data)
    return X, KNN

# Usage exampl


import numpy as np
import pandas as pd
from immupute_dis import immupute_dis
def tfidf(ATAC):
    O = 1 * (ATAC > 0)
    tf1 = O / (np.ones((O.shape[0], 1)) * np.log(1 + np.sum(O, axis=0))[np.newaxis,:])
    idf = np.log(1 + O.shape[1] / (1 + np.sum(O > 0, axis=1)))
    O1 = tf1 * (idf[:, np.newaxis] * np.ones((1, O.shape[1])))
    O1[np.isnan(O1)] = 0
    RE = O1.T
    return RE
def pseudo_bulk(RNA_file,ATAC_file,label_file,Input_dir):
    RNA=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
    Symbol=RNA.index
    RNA = np.log2(1 + RNA.values)
    from sklearn.decomposition import PCA
    # Create an instance of PCA with the desired number of components
    pca = PCA(n_components=100)
    data_pca = pca.fit_transform(RNA.T)
    k = int(np.floor(np.sqrt(RNA.shape[1])))
    TG_filter1, KNN=immupute_dis(RNA.T,data_pca,k)
    ATAC=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    Region=ATAC.index
    ATAC=ATAC.values
    RE=tfidf(ATAC)
    data_pca = pca.fit_transform(RE)
    RE_filter1, KNN=immupute_dis(RE,data_pca,k)
    label=pd.read_csv(Input_dir+label_file,header=None)
    label.columns=['label']
    label=label['label'].values
    cluster=list(set(label))
    A = np.arange(len(label))
    indexfinal=np.array([])
    index0=0
    clusternew=np.array([])
    np.random.seed(42)  # Set seed for reproducibility
    for i in range(len(cluster)):
        d = (label == i)
        idx = np.random.permutation(np.sum(d))
        k = int(np.floor(np.sqrt(np.sum(d))) + 1)
    #print(np.sum(d),k)
    #temp = TG_filter1[:, d]
    #TG_psedubulk[:, index0:(index0+k)] = temp[:, idx[:k]]
        index0 += k
        temp = A[d]
        indexfinal=np.append(indexfinal,temp[:k])
        repeated = np.repeat(cluster[i], k)  # Repeat 0 ten times
        clusternew =np.append(clusternew,repeated)
    indexfinal = indexfinal.astype(int)
    TG_psedubulk = TG_filter1[indexfinal,:].T
    TG_psedubulk=pd.DataFrame(TG_psedubulk,index=Symbol)
    RE_psedubulk = RE_filter1[indexfinal,:].T
    RE_psedubulk=pd.DataFrame(RE_psedubulk,index=Region)
    return TG_psedubulk,RE_psedubulk,indexfinal,clusternew
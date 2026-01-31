import numpy as np
import pandas as pd
import random
import scanpy as sc
#from LingerGRN.immupute_dis import immupute_dis
def tfidf(ATAC):
    O = 1 * (ATAC > 0)
    tf1 = O / (np.ones((O.shape[0], 1)) * np.log(1 + np.sum(O, axis=0))[np.newaxis,:])
    idf = np.log(1 + O.shape[1] / (1 + np.sum(O > 0, axis=1)))
    O1 = tf1 * (idf[:, np.newaxis] * np.ones((1, O.shape[1])))
    O1[np.isnan(O1)] = 0
    RE = O1.T
    return RE
def find_neighbors(adata_RNA,adata_ATAC):
    import scanpy as sc
    K = 20
    #sc.tl.pca(adata_RNA, svd_solver="arpack")
    sc.pp.normalize_total(adata_RNA, target_sum=1e4)
    sc.pp.log1p(adata_RNA)
    sc.pp.highly_variable_genes(adata_RNA, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_RNA.raw=adata_RNA
    adata_RNA = adata_RNA[:, adata_RNA.var.highly_variable]
    sc.pp.scale(adata_RNA, max_value=10)
    sc.tl.pca(adata_RNA, n_comps=15,svd_solver="arpack")
    pca_RNA=adata_RNA.obsm['X_pca']
    sc.pp.log1p(adata_ATAC)
    sc.pp.highly_variable_genes(adata_ATAC, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_ATAC.raw=adata_ATAC
    adata_ATAC = adata_ATAC[:, adata_ATAC.var.highly_variable]
    sc.pp.scale(adata_ATAC, max_value=10, zero_center=True)
    sc.tl.pca(adata_ATAC, n_comps=15,svd_solver="arpack")
    pca_ATAC=adata_ATAC.obsm['X_pca']
    pca = np.concatenate((pca_RNA,pca_ATAC), axis=1)
    adata_RNA.obsm['pca']=pca
    adata_ATAC.obsm['pca']=pca
    #sc.pp.neighbors(adata_RNA, n_neighbors=K, n_pcs=30,use_rep='pca')
    return adata_RNA,adata_ATAC


def pseudo_bulk(adata_RNA,adata_ATAC,singlepseudobulk):
    K = 20
    #sc.tl.pca(adata_RNA, svd_solver="arpack")
    sc.pp.normalize_total(adata_RNA, target_sum=1e4)
    sc.pp.log1p(adata_RNA)
    sc.pp.filter_genes(adata_RNA, min_cells=3) 
    sc.pp.highly_variable_genes(adata_RNA, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_RNA.raw=adata_RNA
    adata_RNA = adata_RNA[:, adata_RNA.var.highly_variable]
    sc.pp.scale(adata_RNA, max_value=10)
    sc.tl.pca(adata_RNA, n_comps=15,svd_solver="arpack")
    pca_RNA=adata_RNA.obsm['X_pca']
    sc.pp.log1p(adata_ATAC)
    sc.pp.filter_genes(adata_ATAC, min_cells=3) 
    sc.pp.highly_variable_genes(adata_ATAC, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_ATAC.raw=adata_ATAC
    adata_ATAC = adata_ATAC[:, adata_ATAC.var.highly_variable]
    sc.pp.scale(adata_ATAC, max_value=10, zero_center=True)
    sc.tl.pca(adata_ATAC, n_comps=15,svd_solver="arpack")
    pca_ATAC=adata_ATAC.obsm['X_pca']
    pca = np.concatenate((pca_RNA,pca_ATAC), axis=1)
    adata_RNA.obsm['pca']=pca
    adata_ATAC.obsm['pca']=pca
    sc.pp.neighbors(adata_RNA, n_neighbors=K, n_pcs=30,use_rep='pca')
    connectivities=(adata_RNA.obsp['distances']>0)
    import random
    label=pd.DataFrame(adata_RNA.obs['label'])
    label.columns=['label']
    label.index=adata_RNA.obs['barcode'].tolist()
    #label=label['label'].values
    cluster=list(set(label['label'].values))
    allindex=[]
    np.random.seed(42)  # Set seed for reproducibility
    for i in range(len(cluster)):
        temp=label[label['label']==cluster[i]].index
        N = len(temp) # Total number of elements
        if N>=10:
            m = int(np.floor(np.sqrt(N)))+1  # Number of elements to sample
            if singlepseudobulk>0:
                m=1
            sampled_elements = random.sample(range(N), m)
            temp=temp[sampled_elements]
            allindex=allindex+temp.tolist()
    connectivities=pd.DataFrame(connectivities.toarray(),index=adata_RNA.obs['barcode'].tolist())
    connectivities=connectivities.loc[allindex].values
    A=(connectivities @ adata_RNA.raw.X.toarray())
    TG_filter1=A/(K-1)
    RE_filter1=(connectivities @ adata_ATAC.raw.X.toarray())/(K-1)
    TG_filter1=pd.DataFrame(TG_filter1.T,columns=allindex,index=adata_RNA.raw.var['gene_ids'].tolist())
    RE_filter1=pd.DataFrame(RE_filter1.T,columns=allindex,index=adata_ATAC.raw.var['gene_ids'].tolist())
    return TG_filter1,RE_filter1

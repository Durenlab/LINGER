def generate_colors(N):
    import matplotlib.colors as mcolors
    import seaborn as sns
    """
    Generate N visually appealing colors using seaborn color palette.
    
    Args:
        N (int): The number of colors to generate.
    
    Returns:
        list: A list of N RGB tuples representing the generated colors.
    """
    color_palette = sns.color_palette("husl", N)
    colors = [mcolors.rgb2hex(color_palette[i]) for i in range(N)]
    return colors
def load_data_ptb(Input_dir,outdir,GRNdir):
    import pandas as pd
    import numpy as np
    import torch
    ATAC_file='ATAC.txt'
    idx_file=outdir+'index.txt'
    RNA_file='RNA.txt'
    label_file='label.txt'
    TFName=outdir+'TFName.txt'
    from LingerGRN import pseudo_bulk
    RNA=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
    ATAC=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    RNA = np.log2(1 + RNA)
    from sklearn.impute import KNNImputer
    K=int(np.floor(np.sqrt(RNA.shape[1])))
    imputer = KNNImputer(n_neighbors=K)
# RNA row is genes col is cells
    TG_filter1 = imputer.fit_transform(RNA.values.T)
    TG_filter1=pd.DataFrame(TG_filter1.T,columns=RNA.columns,index=RNA.index)
    RE_filter1 = imputer.fit_transform(np.log2(1+ATAC.values.T))
    RE_filter1=pd.DataFrame(RE_filter1.T,columns=ATAC.columns,index=ATAC.index)
    #Opn=pd.read_csv(Opn_file,header=0,sep='\t',index_col=0)
    #Opn=Opn.values
    idx=pd.read_csv(idx_file,header=None,sep='\t')
    #Target=pd.read_csv(geneexp_file,header=0,sep='\t',index_col=0)
    genename=pd.read_csv(outdir+'Symbol.txt',sep='\t',header=None)
    genename=genename[0].values
    TFname=pd.read_csv(outdir+'TFName.txt',sep='\t',header=None)
    TFname=TFname[0].values
    Exp=TG_filter1.loc[TFname].values
    Target=TG_filter1.loc[genename].values
    Opn=RE_filter1.values
    chrall=[str(i+1) for i in range(22)]
    chrall.append('X')
    data_merge=pd.read_csv(outdir+'data_merge.txt',sep='\t',index_col=0)
    return chrall,data_merge,Exp,Opn,Target,idx,TFname
def LINGER_simulation(ii,gene_chr,TFindex,Exp,REindex,Opn,netall,index_all):
    import warnings
    import time
    import LingerGRN
    from tqdm import tqdm
    import torch
    import pandas as pd
    import numpy as np
    eps=1e-6
    gene_idx=gene_chr['id_s'].values[ii]-1
    TFidxtemp=TFindex[gene_idx]
    TFidxtemp=TFidxtemp.split('_')
    TFidxtemp=[int(TFidxtemp[k])+1 for k in range(len(TFidxtemp))]
    TFtemp=Exp[np.array(TFidxtemp)-1,:]
    REidxtemp=REindex[gene_idx]
    REidxtemp=str(REidxtemp).split('_')
    if (len(REidxtemp)==1)&(REidxtemp[0]=='nan'):
        REidxtemp=[]
        inputs=TFtemp+1-1
    else:
        REidxtemp=[int(REidxtemp[k])+1 for k in range(len(REidxtemp))]
        REtemp=Opn[np.array(REidxtemp)-1,:]
        inputs=np.vstack((TFtemp, REtemp))
    inputs = torch.tensor(inputs,dtype=torch.float32)
    mean = inputs.mean(dim=1)
    std = inputs.std(dim=1)
    inputs = (inputs.T - mean) / (std+eps)
    inputs=inputs.T              
    num_nodes=inputs.shape[0]   
    loaded_net = netall[index_all]
    X_tr = inputs.T
    y_pred = loaded_net(X_tr)
    return y_pred
def get_simulation(outdir,chrall,data_merge,GRNdir,Exp,Opn,Target,idx):
    import warnings
    import time
    import LingerGRN
    from tqdm import tqdm
    import torch
    import pandas as pd
    import numpy as np
    output=np.zeros(Target.shape)
    for i in range(23):
        chr='chr'+chrall[i]
        print(chr)
        gene_chr=data_merge[data_merge['chr']==chr]
        N=len(gene_chr)
        netall=torch.load(outdir+'net_'+chr+'.pt')
        idx_file1=GRNdir+chr+'_index.txt'
        idx_file_all=GRNdir+chr+'_index_all.txt'
        idxRE_all=pd.read_csv(idx_file_all,header=None,sep='\t')
        gene_chr=data_merge[data_merge['chr']==chr]
        N=len(gene_chr)
        TFindex=idx.values[:,2]
        REindex=idx.values[:,1]
        for ii in tqdm(range(N)):
            index_all=gene_chr.index[ii]
            if index_all in netall.keys():
                res=LINGER_simulation(ii,gene_chr,TFindex,Exp,REindex,Opn,netall,index_all)
                output[index_all,:]=res.detach().numpy().reshape(-1,)
    output1=pd.DataFrame(output,index=data_merge.loc[range(Target.shape[0])]['Symbol'])
    return output1
def umap_embedding(outdir,Target,original,perturb,Input_dir):
    import umap
    import scanpy as sc
# Assuming you have loaded or created an AnnData object named 'adata'
# Create and train the UMAP model
    from sklearn.decomposition import PCA
    import numpy as np
    import pandas as pd
#RNA=pd.read_csv(Input_dir+'RNA.txt',header=0,index_col=0,sep='\t')
    Symbol=pd.read_csv(outdir+'Symbol.txt',header=None,sep='\t')
#sampleall=RNA.columns
    RNA=pd.DataFrame(Target,index=Symbol[0].values)
# Assuming your feature * sample matrix is stored in a variable "matrix"
# Step 1: Calculate the variance across the samples for each feature
    variance = np.var(RNA.values, axis=1)
# Step 2: Sort the features based on the variance in descending order
    sorted_indices = np.argsort(variance)[::-1]
# Step 3: Select the top 2000 features
    top_2000_features = sorted_indices[:2000]
# Assuming you have loaded or created an AnnData object named 'adata'
# Perform PCA using the 'arpack' solver
    pca = PCA(svd_solver='arpack')
    pca.fit(RNA.values[top_2000_features,:].T)
    pca_result = pca.fit_transform(RNA.values[top_2000_features,:].T)
    pca_result=pca_result[:,1:20]
    original1=original.loc[RNA.index[top_2000_features]].values
    original1[original1<0]=0
    original1=np.log2(1+original1)
    #original1=original1
    perturb1=perturb.loc[RNA.index[top_2000_features]].values
    perturb1[perturb1<0]=0
    perturb1=np.log2(perturb1+1)
    #perturb1=perturb1
    O_PCA=pca.fit_transform(original1.T)
    P_PCA=pca.fit_transform(perturb1.T)
    O_PCA=O_PCA[:,1:20]
    P_PCA=P_PCA[:,1:20]
    umap_model = umap.UMAP(n_components=2)
    umap_model.fit(pca_result)
    O_cell_umap = umap_model.transform(O_PCA)
    P_cell_umap = umap_model.transform(P_PCA)
    embedding = umap_model.transform(pca_result)
    D=P_cell_umap-O_cell_umap
    return embedding,D

# Assuming you have a continuous value stored in `continuous_values`
# Define the color for small values (e.g., white) and the color for higher values
def diff_umap(TFko,TFName,save,outdir,embedding,perturb,original,Input_dir):
    import seaborn as sns
    import matplotlib.colors as mcolors
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    label=pd.read_csv(Input_dir+'label.txt',sep='\t',header=None)
    label=label[0].values
    from matplotlib.colors import LinearSegmentedColormap
    zero_color = 'white'
    positive_color = 'orange'
    negative_color = 'blue'
# Define the colormap with a white-to-orange-to-blue gradient
    cmap_colors = [negative_color, zero_color, positive_color]
# Define the colors for each cluster
    sns.set(style='white')
    fig, ax = plt.subplots(figsize=(4, 4))
    continuous_values=perturb.loc[TFName].values-original.loc[TFName].values
#continuous_values[continuous_values<0]=0
    cmap = LinearSegmentedColormap.from_list('custom_cmap', cmap_colors)
# Create a scatter plot with colored dots based on the cluster annotations
    plt.scatter(embedding[:,0], embedding[:,1], c=continuous_values,cmap=cmap,s=2,
            vmin=-np.abs(continuous_values).max(), vmax=np.abs(continuous_values).max())
    anno=label
    unique_clusters = np.unique(anno)
    for cluster in unique_clusters:
        indices = np.where(anno == cluster)
        cluster_center = (np.mean(embedding[indices, 0]), np.mean(embedding[indices, 1]))
        plt.text(cluster_center[0], cluster_center[1], f'Cluster {cluster}', fontsize=10, ha='center', va='center')
    plt.colorbar() 
    plt.xlabel('Umap 1')
    plt.ylabel('Umap 2')
    if save==True:
        plt.savefig(outdir+TFko+"_KO_Diff_exp_Umap_"+TFName+".png", format='png', bbox_inches='tight')
# Add arrows to indicate gene expression changes
    plt.show()
    plt.close()
# Define the colors for each cluster
def Umap_direct(TFko,Input_dir,embedding,D,save,outdir):
    import seaborn as sns
    import matplotlib.colors as mcolors
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt 
    label=pd.read_csv(Input_dir+'label.txt',sep='\t',header=None)
    label=label[0].values
    N=len(np.unique(label))
    colors = generate_colors(N)
    sns.set(style='white')
    label1=label.copy()
    anno=label
    unique_clusters = np.unique(anno)
    D[np.abs(D[:,0])<np.abs(D).mean(axis=0)[0],0]=0
    D[np.abs(D[:,1])<np.abs(D).mean(axis=0)[1],1]=0
    idx=((np.abs(D)>0).sum(axis=1)>=1)
    if type(label[0]) is str:
        for i in range(N):
            label1[label==unique_clusters[i]]=i
    fig, ax = plt.subplots(figsize=(4, 4))
    continuous_values=[colors[i] for i in label1]
# Create a scatter plot with colored dots based on the cluster annotations
    plt.scatter(embedding[:,0], embedding[:,1], c=continuous_values, s=2)
    for cluster in unique_clusters:
        indices = np.where(anno == cluster)
        cluster_center = (np.mean(embedding[indices, 0]), np.mean(embedding[indices, 1]))
        plt.text(cluster_center[0], cluster_center[1], f'Cluster {cluster}', fontsize=10, ha='center', va='center')
# Add arrows to indicate gene expression changes
    ax.quiver(embedding[idx,0], embedding[idx, 1],
           2*D[idx,0],  # Assuming gene_index is the index of the gene you are interested in
          2*D[idx,1],  # Assuming gene_index+1 is the index of another gene for the y-component
          scale=30, scale_units='inches', alpha=0.5)
    if save==True:
        plt.savefig(outdir+TFko+"_KO_Differentiation_Umap.png", format='png', bbox_inches='tight')
    plt.show()
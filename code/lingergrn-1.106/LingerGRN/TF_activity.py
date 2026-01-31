def quantile_normalize(df):
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
def bulk_reg(outdir,GRNdir,genome,chrN):    
    from LingerGRN.LL_net import load_region
    from LingerGRN.LL_net import load_TFbinding
    from LingerGRN.LL_net import load_TF_RE
    from LingerGRN.LL_net import load_RE_TG
    from LingerGRN.LL_net import load_RE_TG_distance
    from scipy.sparse import coo_matrix
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(outdir,GRNdir,genome,chrN)
    TFbinding=load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN)
    mat=load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u)
    TFoverlap = list(set(mat.columns) & set(TFbinding.columns))
    mat = mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    S = np.log(mat+TFbinding+0.1)
    S.index=N_overlap
    S = S.groupby(S.index).max()
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    sparse_S+=0.1
    sparse_dis=load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGset)
    Score=sparse_S.multiply(sparse_dis.values)
    Score=pd.DataFrame(Score.values,index=N_overlap,columns=TGset)
    Score=Score.groupby(Score.index).max()
    data = Score.values[Score.values!=0] 
    rows, cols = np.nonzero(Score.values) 
    coo = coo_matrix((data,(rows,cols)),shape=Score.shape)
    combined = np.zeros([len(data),3], dtype=object) 
    combined[:,0]=Score.index[coo.row]
    combined[:,1]=np.array(TGset)[coo.col]
    combined[:,2]=coo.data
    combined=pd.DataFrame(combined)
    return S,combined
def TF_RE2m(result_RE_TG,REset):
    from scipy.sparse import coo_matrix
    TGset=result_RE_TG['TG'].unique()
    #REset=result_RE_TG.['TG'].unique()
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
# Map the column names and row names to integer indices in the DataFrame
    result_RE_TG["col_index"] = result_RE_TG["TG"].map(col_dict)
    result_RE_TG["row_index"] = result_RE_TG["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = result_RE_TG["col_index"].tolist()
    row_indices = result_RE_TG["row_index"].tolist()
    values = result_RE_TG["Score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)),shape=(len(REset), len(TGset)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
    cis=sparse_S.toarray()
    cis=pd.DataFrame(cis,index=REset,columns=TGset)
    return cis
import scipy.io as sio
import numpy as np
import pandas as pd
def regulon(outdir,adata_RNA,GRNdir,network,genome):
# Load data from MATLAB .mat files
    if network=='cell population':
        trans_reg = pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)
# Apply quantile normalization to 'trans_reg_n'
    elif network=='general':    
        from tqdm import tqdm
        chrom=['chr'+str(i+1) for i in range(22)]
        chrom.append('chrX')
        result_TF_RE=pd.DataFrame([])
        result_RE_TG=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            TF_RE,RE_TG=bulk_reg(outdir,GRNdir,genome,chrN)
            RE_TG.columns=['RE','TG','Score']
            result_RE_TG=pd.concat([result_RE_TG,RE_TG],axis=0,join='outer')
            result_TF_RE=pd.concat([result_TF_RE,TF_RE],axis=0,join='outer')
        TFset=result_TF_RE.columns
        REset=result_TF_RE.index
        cis=TF_RE2m(result_RE_TG,REset)
        TGset=cis.columns
        from LL_net import load_TF_TG
        TF_TG=load_TF_TG(GRNdir, TFset,TGset)
        trans_reg=np.matmul(result_TF_RE.values.T, cis.values).T*(TF_TG.values)
        trans_reg=pd.DataFrame(trans_reg, index=TGset,columns=TFset)
    else: 
        trans_reg = pd.read_csv(outdir+'cell_type_specific_trans_regulatory_'+network+'.txt',sep='\t',index_col=0)
    #RNA=pd.read_csv(outdir+RNA_file,sep='\t',index_col=0)
    RNA=pd.DataFrame(adata_RNA.X.toarray().T,index=adata_RNA.var['gene_ids'].values,columns=adata_RNA.obs['barcode'].values)
    gene_overlap=list(set(trans_reg.index)&set(RNA.index))  
    RNA=RNA.loc[gene_overlap]
    trans_reg=trans_reg.loc[gene_overlap]
    row=trans_reg.sum(axis=1).values[:,np.newaxis]+trans_reg.sum(axis=1).mean()*0.000001
    colsum=trans_reg.sum(axis=0).values[:,np.newaxis]+trans_reg.sum(axis=0).mean()*0.000001
    data_norm=trans_reg.values/(row*colsum.T)*row.sum()
    E=(row*colsum.T)/row.sum()
    
    #trans_reg=trans_reg/trans_reg.sum(axis=0) 
    #trans_reg_norm = quantile_normalize(trans_reg.T)
    RNA=RNA/RNA.sum(axis=0) #neccessary
    RNA_norm=quantile_normalize(RNA)
    #RNA_norm=RNA
    regulon=(np.dot(trans_reg.values.T,RNA_norm.values)-np.dot(E.T,RNA_norm.values))/np.dot(E.T,RNA_norm.values)
    regulon=pd.DataFrame(regulon,index=trans_reg.columns,columns=RNA_norm.columns)
    return regulon
def master_regulator(regulon_score,adata_RNA,celltype):
    import pandas as pd
    import statsmodels.stats.multitest as smm
    from scipy import stats
    #label=pd.read_csv(outdir+labels,sep='\t',header=None)
    #label = label.astype(str)
    #label.columns=['celltype']
    label=pd.DataFrame(adata_RNA.obs['label'].values.tolist(),columns=['celltype'])
    label = label.astype(str)
    if celltype in label['celltype'].values:
# Assuming X and Y are DataFrames with multiple variables/columns
        idx=regulon_score.columns[np.isin(label['celltype'],celltype)]
        X=regulon_score[idx]
        idx=regulon_score.columns[np.isin(label['celltype'],celltype)==0]
        Y=regulon_score[idx]
# Initialize an empty DataFrame to store the t-test results
        t_test_results = np.zeros((regulon_score.shape[0],2))
    # Iterate over the columns/variables in X and Y
        for i in range(X.shape[0]):
            row= regulon_score.index[i]
    # Perform the t-test between X and Y for the current variable
            t_stat, p_value = stats.ttest_ind(X.loc[row], Y.loc[row],alternative='greater')   
    # Append the results to the DataFrame
            t_test_results[i,0] = t_stat 
            t_test_results[i,1] = p_value
        t_test_results[np.isnan(t_test_results)]=1
        t_test_results=pd.DataFrame(t_test_results,index=regulon_score.index,columns=['t_stat','p_value'])
        t_test_results['adj_p']= smm.multipletests(t_test_results['p_value'], method='fdr_bh')[1]
    elif celltype=='all':
        label_set=np.array(list(set(label['celltype'].values)))
        t_test_results = np.zeros((regulon_score.shape[0],3*len(label_set)))
        for j in range(len(label_set)):
            idx=regulon_score.columns[np.isin(label['celltype'],label_set[j])]
            X=regulon_score[idx]
            idx=regulon_score.columns[np.isin(label['celltype'],label_set[j])==0]
            Y=regulon_score[idx]
# Initialize an empty DataFrame to store the t-test results
    # Iterate over the columns/variables in X and Y
            for i in range(X.shape[0]):
                row= regulon_score.index[i]
    # Perform the t-test between X and Y for the current variable
                t_stat, p_value = stats.ttest_ind(X.loc[row], Y.loc[row],alternative='greater')   
    # Append the results to the DataFrame
                p_value = np.nan_to_num(p_value, nan=1)
                t_test_results[i,3*j+1] = p_value 
                t_test_results[i,3*j] = t_stat
            t_test_results[:,3*j+2]=smm.multipletests(t_test_results[:,3*j+1], method='fdr_bh')[1]
        t_test_results=pd.DataFrame(t_test_results,index=regulon_score.index)
        col=[0 for kk in range(len(label_set)*3)]
        for j in range(len(label_set)):
            col[3*j]=label_set[j]+'_t_stat'
            col[3*j+1]=label_set[j]+'_p_value'
            col[3*j+2]=label_set[j]+'_adj_p'
        t_test_results.columns=col
    return t_test_results
def heatmap_cluster(regulon_score,adata_RNA,save,outdir):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import zscore
# Generate random data for the heatmap
    Vars=regulon_score.var(axis=1)
    regulon_score=regulon_score.loc[regulon_score.index[Vars>0]]
    z_scores = zscore(regulon_score, axis=1)
    z_scores1=z_scores.values
    label=adata_RNA.obs['label'].values
    labelset=set(label)
    idx=0
    for labeltemp in labelset:
        index=(label==labeltemp)
        z_scores1[:,idx:idx+index.sum()]=z_scores[z_scores.columns[index]].values
        idx=idx+index.sum()
    z_scores1[z_scores1<-2]=-2
    z_scores1[z_scores1>2]=2
# Set up the figure size
    plt.figure(figsize=(8, 6))
    sns.clustermap(z_scores1, row_cluster=True, col_cluster=False)
    plt.xlabel('Columns')
    plt.ylabel('Rows')
    if save==True:
        plt.savefig(outdir+"heatmap_activity.png", format='png', bbox_inches='tight')
    # Finally, display the plot
    plt.show()
    
def box_comp(TFName,adata_RNA,celltype1,celltype2,datatype,regulon_score,save,outdir):
    import numpy as np
    data=np.zeros(regulon_score.shape[1])
    if datatype=='activity':
        TFexp=regulon_score.loc[TFName].values
    if datatype=='expression':
        data0=pd.DataFrame(adata_RNA.X.toarray().T,index=adata_RNA.var['gene_ids'].values,columns=adata_RNA.obs['barcode'].values)
        TFexp=data0.loc[TFName].values
    label=adata_RNA.obs['label'].values
    if type(label[0]) in [np.int64,np.int32,np.float64,np.float32]:
        label=[str(label[i]) for i in range(len(label))]
    if celltype1=='Others':
        G2=TFexp[np.array(label)==celltype2]
        G1=TFexp[np.array(label)!=celltype2]
    elif celltype2=='Others':
        G1=TFexp[np.array(label)==celltype1]
        G2=TFexp[np.array(label)!=celltype1]
    else:
        G1=TFexp[np.array(label)==celltype1]
        G2=TFexp[np.array(label)==celltype2]
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    # Combine the vectors into a single list
    data = [G1, G2]
    # Set up the figure size and style
    plt.figure(figsize=(4, 3))
    sns.set_style("white")
    # Create a list of colors for each violin
    colors = ["#b2103e", "#22b6ed"]
    # Plot the violin plot with different colors for each violin
    sns.violinplot(data=data, palette=colors)
    # Add box plots to the violin plot
    sns.boxplot(data=data, color="white", width=0.15)
# Customize the plot
    plt.xlabel('Violins')
    plt.ylabel('Values')
# Rename x-axis labels
    plt.xticks([0, 1], [celltype1, celltype2])
    if save==True:
        plt.savefig(outdir+"box_plot_"+TFName+'_'+datatype+'_'+celltype1+'_'+celltype2+".png", format='png', bbox_inches='tight')
# Finally, display the plot
    plt.show()
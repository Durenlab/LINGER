from scipy import stats
import pandas as pd
import numpy as np
def assignLabel(W,p):
    W=W/(W.sum(axis=0)+1**(-6))
    W2=W.T/(W.T.sum(axis=0)+ 10**(-4));
    max_values = np.max(W2, axis=0)#col max
    max_indices = np.argmax(W2, axis=0)
    quantile = np.percentile(max_values, p*100)
    S_gene=W[:,0];
    S_gene[:]=0
    K=W.shape[1]
    for i in range(K):
        S_gene[(max_values>quantile)&(max_indices==i)]=i+1
    return S_gene,W2 
import numpy as np
import matplotlib.pyplot as plt

def qq_pval(p1,names,celltype):
    p1[np.isnan(p1)] = 1
    d = np.sort(p1)
    f = np.argsort(p1)
    p11 = np.arange(1, len(p1) + 1) / len(p1)
    plt.scatter(-np.log10(p11), -np.log10(d))
    for i in range(K):
        plt.text(-np.log10(p11[i])+np.log10(K)/20, -np.log10(d[i]), np.array(names)[f[i]], fontsize=8)

    m = round(-np.log10(d[0]) * 100) / 100
    plt.plot(np.arange(0, m + 0.01, 0.01), np.arange(0, m + 0.01, 0.01), '-r', linewidth=1)
    plt.xlabel('Theoretical Quantiles')
    plt.ylabel('Sample Quantiles')
    plt.title(celltype)
    plt.show()  # Display the plot
    plt.close()

class Module_obj:
    def __init__(self):
        import pandas as pd
        import numpy as np
        self.S_TG = pd.DataFrame()  # Initialize A.x as an empty DataFrame
        self.pvalue_all = pd.DataFrame()  # Initialize A.y as an empty list
        self.tvalue_all = pd.DataFrame() 
        self.p_fisher = pd.DataFrame() 
        self.odds_ratio_fisher=pd.DataFrame() 
        

def diff_Module(Exp_TG,metadata,S_TG,K):
    celltype=metadata['celltype'].unique().tolist()
    pvalue_all= np.zeros((K, len(celltype)))
    tvalue_all= np.zeros((K, len(celltype)))
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    for k in range(len(celltype)):
        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype[k]]
        aud_idxtemp=metadata[(metadata['celltype'].values==celltype[k])]['group'].values
        Exp_mean=stats.zscore(temp.T).T.groupby(S_TG['Module'].values).mean()
        Exp_mean=Exp_mean.loc[range(1,K+1)]
        X=Exp_mean.values[:,(aud_idxtemp==1)]
        Y=Exp_mean.values[:,(aud_idxtemp==0)]
        p_values = np.zeros((K, ))   
        t_values = np.zeros((K, ))
        from scipy.stats import ttest_ind
        for i in range(K):
            t_values[i], p_values[i] = ttest_ind(X[i], Y[i])
        #p_values = np.nan_to_num(p_values, nan=1)
        pvalue_all[:,k]=p_values
        tvalue_all[:,k]=t_values
    
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
    pvalue_all=pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype)
    tvalue_all=pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype)
    return pvalue_all,tvalue_all


def GWAS_Module_enrich(S_TG,TGset,GWASgene,K):
    from collections import Counter
    counts = Counter(S_TG['Module'].values)
    from scipy.stats import fisher_exact
    p_fisher=np.zeros((K,GWASgene.shape[1]));
    odds_ratio_fisher=np.zeros((K,GWASgene.shape[1]));
    for i in range(GWASgene.shape[1]):
        AUD_genes=GWASgene[GWASgene['GWAS_'+(str(i+1))]==1].index
        N=len(set(TGset)&set(AUD_genes))
        for k in range(K):
            Noverlap=len(set(TGset[S_TG['Module'].values==k+1])&set(AUD_genes))
            contingency_table = np.array([[Noverlap,  counts[k+1]-Noverlap], [N-Noverlap, len(TGset)- N-counts[k+1]+Noverlap]])
            odds_ratio, p_value = fisher_exact(contingency_table,alternative='greater')
            p_fisher[k,i]=odds_ratio
            odds_ratio_fisher[k,i]=p_value
    odds_ratio_fisher=pd.DataFrame(odds_ratio_fisher,index=['M'+str(j+1) for j in range(K)],columns=['GWAS_'+(str(j+1)) for j in range(GWASgene.shape[1])])
    p_fisher=pd.DataFrame(p_fisher,index=['M'+str(j+1) for j in range(K)],columns=['GWAS_'+(str(j+1)) for j in range(GWASgene.shape[1])])
    return p_fisher,odds_ratio_fisher

def Module_trans(outdir,metadata,TG_pseudobulk,K,GWASfile=None):
    import numpy as np
    from scipy import stats
    print('loading GRN......')
    trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)
    #print(trans_reg)
    #TG_pseudobulk=TG_pseudobulk_all
    TFset=trans_reg.columns
    TGset=trans_reg.index
    #TG_pseudobulk=TG_pseudobulk/TG_pseudobulk.mean()
    idx = [s[:3] != "MIR" for s in TG_pseudobulk.index]
    TG_pseudobulk=TG_pseudobulk.loc[TG_pseudobulk.index[idx]]
    R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
    R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
    #R1[R1<0]=0;R2[R2<0]=0
    #X=TG_pseudobulk.loc[TGset]
    from sklearn.preprocessing import quantile_transform
    #Exp=quantile_transform(np.log2(X + 1), n_quantiles=10, random_state=0, copy=True)
    #Exp=quantile_transform(np.log2(TG_pseudobulk + 1), n_quantiles=10, random_state=0, copy=True)
    Exp=pd.DataFrame(TG_pseudobulk,index=TG_pseudobulk.index,columns=TG_pseudobulk.columns)
    Z=R1+R2
    Z[Z<0]=0  
    print('identify modules......')
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=K, init='random', random_state=0)
    W = nmf.fit_transform(Z)
    H = nmf.components_
    [S_TG,W2]=assignLabel(W,0.9);
    [S_TF,H2]=assignLabel(H.T,0.8);
    #Exp_mean=stats.zscore(Exp_TG.T).T.groupby(S_TG).mean()
    S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
    Exp_TF=Exp.loc[TFset]
    Exp_TG=Exp.loc[TGset]
    print('differential modules......')
    pvalue_all,tvalue_all=diff_Module(Exp_TG,metadata,S_TG,K)
    if GWASfile is not None:            
        print('GWAS enrich......')
        GWASgene=pd.DataFrame(index=TG_pseudobulk.index)
        for i in range(len(GWASfile)):
            temp=pd.read_csv(GWASfile[i],sep='\t',header=None)
            idx=np.zeros((GWASgene.shape[0],1))
            idx[GWASgene.index.isin(temp[0].values),:]=1
            GWASgene['GWAS_'+(str(i+1))]=idx
        p_fisher,odds_ratio_fisher=GWAS_Module_enrich(S_TG,TGset,GWASgene,K)
        Module_result=Module_obj()
        Module_result.S_TG=S_TG
        Module_result.pvalue_all=pvalue_all
        Module_result.tvalue_all=tvalue_all
        Module_result.p_fisher=p_fisher
        Module_result.odds_ratio_fisher=odds_ratio_fisher
        return Module_result
    else:
        Module_result=Module_obj()
        Module_result.S_TG=S_TG
        Module_result.pvalue_all=pvalue_all
        Module_result.tvalue_all=tvalue_all
        return Module_result

def remove_covariate(TG_pseudobulk,aud_idx,celltype):
    celltypeset=list(set(celltype))
    Exp_norm=TG_pseudobulk.values.copy()
    for i in range(len(celltypeset)):
        Nseq=np.array(range(len(celltype)))
        Cov=aud_idx.iloc[Nseq[celltype==celltypeset[i]],:]
        Cov_product = np.dot(Cov.T, Cov)
        # Calculate the pseudo-inverse of Cov_product
        Cov_product_pinv = np.linalg.pinv(Cov_product)
        # Calculate beta
        Exp=TG_pseudobulk.iloc[:,Nseq[celltype==celltypeset[i]]]
        beta = np.dot(np.dot(Exp, Cov), Cov_product_pinv)
        Exp_norm1=Exp-np.dot(beta, Cov.T)+np.dot(np.reshape(beta[:, 0], (-1, 1)),np.ones((1,Exp.shape[1])))
        Exp_norm1[Exp_norm1<0]=0
        Exp_norm[:,celltype==celltypeset[i]]=Exp_norm1
    Exp_norm=pd.DataFrame(Exp_norm,columns=TG_pseudobulk.columns,index=TG_pseudobulk.index)
    return Exp_norm

import pandas as pd
import numpy  as np
def runWGCNA(celltypetemp,TG_pseudobulk_all,metadata):
    metadata_temp=metadata[metadata['celltype'].isin([celltypetemp])]
    idx=np.array([i for i in range(metadata.shape[0])])
    idx=idx[metadata['celltype'].isin([celltypetemp])]
    TG_pseudobulk=TG_pseudobulk_all.iloc[:,idx]
    print(TG_pseudobulk.shape)
    TG_pseudobulk=TG_pseudobulk.groupby(TG_pseudobulk.index).mean()
    print(TG_pseudobulk.shape)
    expression=TG_pseudobulk.T
    print(expression.shape)
    geneList=pd.DataFrame(TG_pseudobulk.index,index=TG_pseudobulk.index)
    geneList.columns=['gene_name']
    geneList['gene_type']='protein_coding'
    import PyWGCNA
    expression.to_csv(celltypetemp+'_expression_WGCNA.csv')
    pyWGCNA_5xFAD = PyWGCNA.WGCNA(name=celltypetemp,  
                              geneExpPath=celltypetemp+'_expression_WGCNA.csv', 
                              outputPath='',
                              save=True)
    pyWGCNA_5xFAD.preprocess()
    pyWGCNA_5xFAD.findModules()
    metadata_temp.to_csv(celltypetemp+'_sampleInfo.csv',sep=',')
    pyWGCNA_5xFAD.updateSampleInfo(path=celltypetemp+'_sampleInfo.csv', sep=',')
    pyWGCNA_5xFAD.setMetadataColor('group', {0: 'green',
                                       1: 'yellow'})
    pyWGCNA_5xFAD.updateGeneInfo(geneList)
    pyWGCNA_5xFAD.analyseWGCNA()
    pyWGCNA_5xFAD.saveWGCNA()

def correlation_FC(x,y,method):
    import numpy as np
    from scipy import stats
    # Loop through each column of y and calculate correlation with x
    correlations = []
    correlationsp=[]
    for i in range(y.shape[1]):
        if method=='pearsonr':
            r, p = stats.pearsonr(x.ravel(), y.values[:,i])  
        if method=='spearmanr':
            r, p = stats.spearmanr(x.ravel(), y.values[:,i]) 
        correlations.append(r)
        correlationsp.append(p)
    correlations=pd.DataFrame(correlations,index=y.columns)
    correlationsp=pd.DataFrame(correlationsp,index=y.columns)
    return correlations,correlationsp


def driver_score(expression,aud_idx,GRN,outdir,adjust_method,corr_method):
    print('loading GRN......')
    import numpy as np   
    from statsmodels.stats.multitest import multipletests
    allcelltype=aud_idx['celltype'].unique()
    C_result=pd.DataFrame([])
    P_result=pd.DataFrame([])
    Q_result=pd.DataFrame([])
    reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)
    #print(reg.shape)
    if reg.shape[1]<4:
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    overlap=list(set(expression.index)&set(reg.index))
    reg=reg.loc[overlap]
    expression=expression.loc[overlap]
    cols=expression.sum(axis=0).values
    rows=expression.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    E=E+E.mean().mean()*10**(-1)
    expression=expression/E
    for i in range(len(allcelltype)):
        print('cell type '+ allcelltype[i])
        aud_idx1=aud_idx.reset_index()
        exp_temp=expression.iloc[:,aud_idx1[aud_idx1['celltype']==allcelltype[i]].index]
        aud_idx1=aud_idx1[aud_idx1['celltype']==allcelltype[i]]
        Mean1=exp_temp[aud_idx1[aud_idx1['group']==1]['index']].mean(axis=1)+10**(-6)
        Mean0=exp_temp[aud_idx1[aud_idx1['group']==0]['index']].mean(axis=1)+10**(-6)
        FC=pd.DataFrame(Mean1/Mean0)
        FC=FC.loc[reg.index]
        #print(FC)
        print(np.isnan(FC).sum())
        c,cp=correlation_FC(np.log(FC[0]).values,reg,corr_method)
        #idx=pd.DataFrame(range(expression.shape[0]),index=expression.index)
        C_result=pd.concat([C_result,c],axis=1)
        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        #print(adjusted_p_values)
        Q_result=pd.concat([Q_result,adjusted_p_values],axis=1)  
    C_result.columns=allcelltype
    P_result.columns=allcelltype
    Q_result.columns=allcelltype    
    return C_result,P_result,Q_result

def driver_result(C_result,Q_result,K):
    import pandas as pd
    # Create a sample DataFrame of size 100x7   
    # Rank all values in the DataFrame
    top_5_rows=[]
    for column in C_result.columns:
        sorted_column = C_result[column].sort_values(ascending=False)
        top_5_rows =top_5_rows+sorted_column.index[:K].tolist()
    top_5_rows=list(set(top_5_rows))
    for column in C_result.columns:
        sorted_column = C_result[column].sort_values(ascending=True)
        top_5_rows =top_5_rows+sorted_column.index[:K].tolist()
    top_5_rows=list(set(top_5_rows))
    return C_result.loc[top_5_rows],Q_result.loc[top_5_rows]


import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix
from tqdm import tqdm
import torch
import csv
import torch.nn as nn
import torch.optim as optim
from torch.nn import functional as F
from scipy.stats import pearsonr
from scipy.stats import spearmanr
#load data
import random
from torch.optim import Adam
import os
from sklearn.linear_model import ElasticNet
from sklearn.datasets import make_regression
from sklearn.model_selection import KFold
hidden_size  = 64
hidden_size2 = 16
output_size = 1
seed_value = 42
torch.manual_seed(seed_value)
class Net(nn.Module):
    def __init__(self,input_size,activef):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.fc2 = nn.Linear(64, 16)
        self.fc3 = nn.Linear(16, output_size)
        self.activef=activef
    def forward(self, x):
        #x = torch.sigmoid(self.fc1(x))
        if self.activef=='ReLU':
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))
        if self.activef=='sigmoid':
            x = F.sigmoid(self.fc1(x))
            x = F.sigmoid(self.fc2(x))
        if self.activef=='tanh':
            x = F.tanh(self.fc1(x))
            x = F.tanh(self.fc2(x))
        x = self.fc3(x)
        return x

def cosine_similarity_0(X):
    A=X.T/((X**2).sum(axis=1)**(1/2)+((X**2).sum(axis=1)**(1/2)).mean()/1000000)
    return np.dot(A.T,A)

    
def list2mat(df,i_n,j_n,x_n):
    TFs = df[j_n].unique()
    REs = df[i_n].unique()
#Initialize matrix as numpy array 
#Map row and col indices for lookup
    row_map = {r:i for i,r in enumerate(REs)}
    col_map = {c:i for i,c in enumerate(TFs)}
    row_indices = np.array([row_map[row] for row in df[i_n]])
    col_indices = np.array([col_map[col] for col in df[j_n]])
    from scipy.sparse import coo_matrix
    matrix = coo_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))
    mat=coo_matrix.toarray(matrix)
    return mat,REs,TFs



def list2mat_s(df,REs,TFs,i_n,j_n,x_n):
#Initialize matrix as numpy array 
#Map row and col indices for lookup
    row_map = {r:i for i,r in enumerate(REs)}
    col_map = {c:i for i,c in enumerate(TFs)}
    row_indices = np.array([row_map[row] for row in df[i_n]])
    col_indices = np.array([col_map[col] for col in df[j_n]])
    from scipy.sparse import coo_matrix
    import scipy.sparse as sp
    matrix = sp.csr_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))
    return matrix,REs,TFs

def merge_columns_in_bed_file(file_path,startcol):
    merged_values = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1+startcol]
            col2 = columns[startcol]
            col3 = columns[1+startcol]
            merged_value = f"{col1}:{col2}-{col3}"
            merged_values.append(merged_value)
    return merged_values
def merge_columns_in_bed_file2(file_path,startcol):
    merged_values = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1+startcol]
            col2 = columns[startcol]
            col3 = columns[1+startcol]
            merged_value = f"{col1}_{col2}_{col3}"
            merged_values.append(merged_value)
    return merged_values
def format_RE_tran12(region):
    chr, range_ = region.split(":")
    start, end = range_.split("-") 
    return "_".join([chr, start, end])
def get_TF_RE(data_merge_temp,j,net_all,TFindex,TFName,REindex,REName):
    index_all=data_merge_temp[j]
    result={'TF':[],'RE':[],'score':[]}
    result = pd.DataFrame(result)
    #for ii in range(1):
    temps=list(net_all[index_all].parameters())[0]
    TFidxtemp=TFindex[index_all]
    TFidxtemp=TFidxtemp.split('_')
    TFidxtemp=[int(TFidxtemp[i]) for i in range(len(TFidxtemp))]
    TFName_temp=TFName[np.array(TFidxtemp)] 
    REidxtemp=REindex[index_all]
    if REidxtemp=='':
        REidxtemp=[]
    else:
        REidxtemp=REidxtemp.split('_')
        REidxtemp=[int(REidxtemp[i]) for i in range(len(REidxtemp))] #146 RE idx   
    if len(REidxtemp)>0:
        corr_matrix = cosine_similarity_0(temps.detach().numpy().T)
        REName_temp=REName[np.array(REidxtemp)]
        corr_matrix=corr_matrix[:len(TFidxtemp),len(TFidxtemp):]
        for k in range(len(REidxtemp)):
            datatemp=pd.DataFrame({'score':corr_matrix[:,k].tolist()})
            datatemp['TF']=TFName_temp.tolist()
            datatemp['RE']=REName_temp[k]
            result=pd.concat([result,datatemp])    
    return result
def load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN):
    TFbinding=pd.read_csv(GRNdir+'TF_binding_'+chrN+'.txt',sep='\t',index_col=0)
    TFbinding1=np.zeros([len(O_overlap_u),TFbinding.shape[1]])
    TFbinding1=np.zeros([len(O_overlap_u),TFbinding.shape[1]])
    O_overlap1=list(set(O_overlap_hg19_u)&set(TFbinding.index))
    List=pd.DataFrame(range(len(TFbinding.index)), index=TFbinding.index)
    index0=List.loc[O_overlap1][0].values
    #O_overlap_df=pd.DataFrame(range(len(O_overlap)), index=O_overlap)
    O_overlap_hg19_u_df=pd.DataFrame(range(len(O_overlap_hg19_u)), index=O_overlap_hg19_u)
    #hg19_38=pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    index1=O_overlap_hg19_u_df.loc[O_overlap1][0].values
    #index1=O_overlap_df.loc[index1][0].values
    TFbinding1[index1,:]=TFbinding.iloc[index0,:].values
    #TFbinding=pd.DataFrame(TFbinding1,index=O_overlap_u,columns=TFbinding.columns)
    O_overlap_u_df=pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38=pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    TFbinding2=np.zeros([len(O_overlap),TFbinding.shape[1]])
    index=O_overlap_u_df.loc[O_overlap][0].values
    TFbinding2=TFbinding1[index,:]
    TFbinding=pd.DataFrame(TFbinding2,index=O_overlap,columns=TFbinding.columns)
    return TFbinding

def load_region(GRNdir,genome,chrN,outdir):  
    O_overlap=merge_columns_in_bed_file(outdir+'Region_overlap_'+chrN+'.bed',1)
    N_overlap=merge_columns_in_bed_file(outdir+'Region_overlap_'+chrN+'.bed',4)
    O_overlap_u=list(set(O_overlap))
    N_overlap_u=list(set(N_overlap))
    #O_all=merge_columns_in_bed_file(GRNdir+'Peaks_'+chrN+'.bed',1)
    hg19_region=merge_columns_in_bed_file(GRNdir+'hg19_Peaks_'+chrN+'.bed',1)
    hg19_region=pd.DataFrame(range(len(hg19_region)),index=hg19_region)
    hg38_region=merge_columns_in_bed_file(GRNdir+'hg38_Peaks_'+chrN+'.bed',1)
    hg38_region=pd.DataFrame(range(len(hg38_region)),index=hg38_region)
    if genome=='hg19':
        idx=hg19_region.loc[O_overlap_u][0].values
        O_overlap_u=hg38_region.index[idx].tolist()
        O_overlap_hg19_u=hg19_region.index[idx].tolist()
    if genome=='hg38':
        idx=hg38_region.loc[O_overlap_u][0].values
        O_overlap_hg19_u=hg19_region.index[idx].tolist()
    return O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u

def load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u):      
    #print('load prior TF-RE for '+chrN+'...')
    mat=pd.read_csv(GRNdir+'Primary_TF_RE_'+chrN+'.txt',sep='\t',index_col=0)
    mat1=np.zeros([len(O_overlap_u),mat.shape[1]])
    O_overlap1=list(set(O_overlap_u)&set(mat.index))
    List=pd.DataFrame(range(len(mat.index)), index=mat.index)
    index0 = List.loc[O_overlap1][0].values
    O_overlap_u_df=pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    index1 = O_overlap_u_df.loc[O_overlap1][0].values
    mat1[index1,:] = mat.iloc[index0,:].values
    #mat = pd.DataFrame(mat1,index=O_overlap_u,columns=mat.columns)
    O_overlap_u_df=pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38=pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    mat2=np.zeros([len(O_overlap),mat.shape[1]])
    index=O_overlap_u_df.loc[O_overlap][0].values
    mat2=mat1[index,:]
    mat=pd.DataFrame(mat2,index=O_overlap,columns=mat.columns)
    return mat
def TF_RE_LINGER_chr(chr,outdir):
    REName = 'data/Peaks.txt'
# Open the file in read mode
    with open(REName, "r") as file:
    # Create a CSV reader
        reader = csv.reader(file, delimiter='\t')
    # Read the first column and store it in a list
        first_column = [row[0] for row in reader]
    REName=np.array(first_column)
    idx_file=outdir+'index.txt'
    from scipy.stats import zscore
    data0=pd.read_csv(outdir+'result_'+chr+'.txt',sep='\t')
    data0.columns=['gene','x','y']
    idx_file=outdir+'index.txt'
    idx=pd.read_csv(idx_file,sep='\t',header=None)
    idx.columns=['gene','REid','TF_id','REid_b']
    idx.fillna('', inplace=True)
    TFName=outdir+'TFName.txt'
    TFName=pd.read_csv(TFName,sep='\t',header=None)
    TFName.columns=['Name']
    TFName=TFName['Name'].values
    TFindex=idx['TF_id'].values
    REindex=idx['REid'].values
    geneName=idx['gene'].values
    net_all=torch.load(outdir+"net_"+chr+".pt")
    data_merge=pd.read_csv(outdir+'data_merge.txt',sep='\t',header=0,index_col=0)
    data_merge_temp=data_merge[data_merge['chr']==chr].index
    batchsize=50
    AAA=np.abs(data0[['x']].values)
    N=data_merge_temp.shape[0]
    times=int(np.floor(N/batchsize))
    resultlist=[0 for i in range(times+1)]
    for ii in tqdm(range(times)):
        result_all=pd.DataFrame([])
        for j in range(ii*batchsize,(ii+1)*batchsize):
            if (AAA[j]>0)&(AAA[j]<10):
                result=get_TF_RE(data_merge_temp,j,net_all,TFindex,TFName,REindex,REName)
                result_all=pd.concat([result_all,result],axis=0)
        result_all=result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
        resultlist[ii]=result_all
    result_all=pd.DataFrame([])
    ii=ii+1
    for j in range(ii*batchsize,N):
        if (AAA[j]>0)&(AAA[j]<10):
            result=get_TF_RE(data_merge_temp,j,net_all,TFindex,TFName,REindex,REName)
            result_all=pd.concat([result_all,result],axis=0)
    if result_all.shape[0]>0:
        result_all=result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
    resultlist[ii]=result_all
    result_all1=pd.concat(resultlist,axis=0)
    A=result_all1.groupby(['TF', 'RE'])['score'].max().reset_index()
    mat,REs,TFs=list2mat(A,'RE','TF','score')
    mat=pd.DataFrame(mat,index=REs,columns=TFs)
    return mat

def TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,outdir):
    ## the regions
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(GRNdir,genome,chrN,outdir)
    import numpy as np
    import pandas as pd
## read the count file.
    #RE=pd.DataFrame(adata_ATAC.raw.X.toarray().T,index=adata_ATAC.raw.var['gene_ids'].values,columns=adata_ATAC.obs['barcode'].values)
    TG=pd.DataFrame(adata_RNA.X.toarray().T,index=adata_RNA.var['gene_ids'].values,columns=adata_RNA.obs['barcode'].values)
## cell annotation
## extact the overlapped peaks.
    #RE=RE.loc[N_overlap]
    #TFbinding=load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN)
    mat=load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u)
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat[TFoverlap]
    #TFbinding = TFbinding[TFoverlap]
    #TF = TG.loc[TFoverlap]
    #mat_m=np.mean(mat.values[mat>0])
    #mat = mat / mat_m
    mat.values[mat.values<0]=0
    #TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    #TF_cluster = TF.values.mean(axis=1)
    #TF_cluster = TF_cluster[None,:]
    #RE_cluster = RE.values.mean(axis=1)
    #RE_cluster = RE_cluster[:,None]
    #S = np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1) + np.log(TF_cluster+0.1)
    S = mat#+TFbinding
    #S = np.exp(S)
    S.index=N_overlap
    mean_S = S.groupby(S.index).max()
    return mean_S
import ast
def TF_RE_scNN(TFName,geneName,net_all,RE_TGlink,REName):
    batchsize=50
    REName=pd.DataFrame(range(len(REName)),index=REName)
    N=RE_TGlink.shape[0]
    times=int(np.floor(N/batchsize))
    resultlist=[0 for i in range(times+1)]
    for ii in range(times):
        result_all=pd.DataFrame([])
        for j in range(ii*batchsize,(ii+1)*batchsize):
            RE_TGlink_temp=RE_TGlink.values[j,:]
            temps=list(net_all[j].parameters())[0]
            actual_list = ast.literal_eval(RE_TGlink_temp[1])
            REidxtemp=REName.loc[actual_list].index
            TFidxtemp=np.array(range(len(TFName)))
            TFidxtemp=TFidxtemp[TFName!=RE_TGlink_temp[0]]
            if len(REidxtemp)>0:
                corr_matrix = cosine_similarity_0(temps.detach().numpy().T)
                corr_matrix=corr_matrix[:len(TFidxtemp),len(TFidxtemp):]
                result={'TF':[],'RE':[],'score':[]}
                result = pd.DataFrame(result)
                for k in range(len(REidxtemp)):
                    datatemp=pd.DataFrame({'score':corr_matrix[:,k].tolist()})
                    datatemp['TF']=TFName[TFidxtemp].tolist()
                    datatemp['RE']=REidxtemp[k]
                    result=pd.concat([result,datatemp]) 
                result_all=pd.concat([result_all,result],axis=0)
        result_all=result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
        #print(result_all)
        resultlist[ii]=result_all
    result_all=pd.DataFrame([])
    ii=times
    if N>ii*batchsize:
        for j in range(ii*batchsize,N):
            RE_TGlink_temp=RE_TGlink.values[j,:]
            temps=list(net_all[j].parameters())[0]
            actual_list = ast.literal_eval(RE_TGlink_temp[1])
            REidxtemp=REName.loc[actual_list].index
            TFidxtemp=np.array(range(len(TFName)))
            TFidxtemp=TFidxtemp[TFName!=RE_TGlink_temp[0]]
            if len(REidxtemp)>0:
                corr_matrix = cosine_similarity_0(temps.detach().numpy().T)
                corr_matrix=corr_matrix[:len(TFidxtemp),len(TFidxtemp):]
                result={'TF':[],'RE':[],'score':[]}
                result = pd.DataFrame(result)
                for k in range(len(REidxtemp)):
                    datatemp=pd.DataFrame({'score':corr_matrix[:,k].tolist()})
                    datatemp['TF']=TFName[TFidxtemp].tolist()
                    datatemp['RE']=REidxtemp[k]
                    result=pd.concat([result,datatemp]) 
                result_all=pd.concat([result_all,result],axis=0)
        result_all=result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
    #print(result_all)
        resultlist[ii]=result_all
        result_all=pd.concat(resultlist,axis=0)
        result_all=result_all.groupby(['TF', 'RE'])['score'].max().reset_index()
    return result_all
    
def load_data_scNN(GRNdir,genome):
    import pandas as pd
    genome_map=pd.read_csv(GRNdir+'genome_map_homer.txt',sep='\t',header=0)
    genome_map.index=genome_map['genome_short'].values 
    if genome in genome_map.index:
        Match2=pd.read_csv(GRNdir+'Match_TF_motif_'+genome_map.loc[genome]['species_ensembl']+'.txt',sep='\t',header=0)  
    else:
        Match2=pd.read_csv(GRNdir+'MotifMatch.txt',sep='\t',header=0)
    TFName = pd.DataFrame(Match2['TF'].unique())
    Target=pd.read_csv('data/TG_pseudobulk.tsv',sep=',',header=0,index_col=0)
    TFlist=list(set(Target.index)&set(TFName[0].values))
    Exp=Target.loc[TFlist]
    Opn=pd.read_csv('data/RE_pseudobulk.tsv',sep=',',header=0,index_col=0)
    RE_TGlink=pd.read_csv('data/RE_gene_distance.txt',sep='\t',header=0)
    RE_TGlink = RE_TGlink.groupby('gene').apply(lambda x: x['RE'].values.tolist()).reset_index()
    geneoverlap=list(set(Target.index)&set(RE_TGlink['gene']))
    RE_TGlink.index=RE_TGlink['gene']
    RE_TGlink=RE_TGlink.loc[geneoverlap]
    RE_TGlink=RE_TGlink.reset_index(drop=True)
    return Exp,Opn,Target,RE_TGlink  
    
    
def TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir):
    from tqdm import tqdm
    import numpy as np
    import pandas as pd
    print('Generating cellular population TF binding strength ...')
    chrom = ['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    if method=='baseline':
        result=pd.DataFrame()
        for i in tqdm(range(23)):
            chrN=chrom[i]
            out=TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,outdir)
            out.to_csv(outdir+chrN+'_cell_population_TF_RE_binding.txt',sep='\t')  
        #result=pd.concat([result,out],axis=1).fillna(0)
            result = pd.concat([result, out], join='outer', axis=0)
    if method=='LINGER':
        result=pd.DataFrame()
        for i in tqdm(range(23)):
            chrN=chrom[i]
            print('Generating cellular population TF binding strength for '+chrN)
            mat=TF_RE_LINGER_chr(chrN,outdir)
            TFs = mat.columns
## read the count file.
            TG=pd.DataFrame(adata_RNA.X.toarray().T,index=adata_RNA.var['gene_ids'].values,columns=adata_RNA.obs['barcode'].values)
            TFoverlap = list(set(TFs) & set(TG.index))
            mat = mat[TFoverlap]
            mat.to_csv(outdir+chrN+'_cell_population_TF_RE_binding.txt',sep='\t')  
            result = pd.concat([result, mat], join='outer', axis=0)
    if method=='scNN':
        Exp,Opn,Target,RE_TGlink=load_data_scNN(GRNdir,genome)
        RE_TGlink=pd.read_csv(outdir+'RE_TGlink.txt',sep='\t',header=0)
        RE_TGlink.columns=[0,1,'chr']
        #chrall=[RE_TGlink[0][i][0].split(':')[0] for i in range(RE_TGlink.shape[0])]
        chrlist=RE_TGlink['chr'].unique()
        REName=Opn.index
        geneName=Target.index
        TFName=Exp.index
        result_all=pd.DataFrame([])
        for jj in tqdm(range(0,len(chrlist))):
            chrtemp=chrlist[jj]
            RE_TGlink1=RE_TGlink[RE_TGlink['chr']==chrtemp]
            net_all=torch.load(outdir+chrtemp+'_net.pt')
            result=TF_RE_scNN(TFName,geneName,net_all,RE_TGlink1,REName)
            result.to_csv(outdir+chrtemp+'_cell_population_TF_RE_binding.txt',sep='\t')
            result_all=pd.concat([result_all,result],axis=0)
        result=result_all.copy()
    result.to_csv(outdir+'cell_population_TF_RE_binding.txt',sep='\t')   
        
def load_TFbinding_scNN(GRNdir,outdir,genome):
    import pandas as pd
    import numpy as np
    genome_map=pd.read_csv(GRNdir+'genome_map_homer.txt',sep='\t',header=0)
    genome_map.index=genome_map['genome_short'].values 
    A=pd.read_csv(outdir+'MotifTarget.bed',sep='\t',header=0,index_col=None)     
            #Motif_binding,REs1,motifs=list2mat(A,'PositionID','Motif Name','MotifScore')
    A['MotifScore']=np.log(1+A['MotifScore']); 
    if genome in genome_map.index:
        Match2=pd.read_csv(GRNdir+'Match_TF_motif_'+genome_map.loc[genome]['species_ensembl']+'.txt',sep='\t',header=0)
    else:
        Match2=pd.read_csv(GRNdir+'MotifMatch.txt',sep='\t',header=0)
    TF_binding,REs1,motifs=list2mat(A,'PositionID','Motif Name','MotifScore')
    TF_binding1=pd.DataFrame(TF_binding.T,index=motifs,columns=REs1)
    TF_binding1['motif']=motifs
    TF_binding1=TF_binding1.merge(Match2,how='inner',left_on='motif',right_on='Motif')
    TF_binding=TF_binding1.groupby(['TF'])[REs1].max()
    TF_binding=TF_binding.reset_index()
    TF_binding.index=TF_binding['TF']
    TF_binding=TF_binding[REs1]
    return TF_binding.T

def cell_type_specific_TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,celltype,outdir,method,mat):
    ## the regions
        ## the regions
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(GRNdir,genome,chrN,outdir)
    import numpy as np
    import pandas as pd
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    temp=adata_ATAC.X[np.array(label)==celltype,:].mean(axis=0)
    RE=pd.DataFrame(temp.T,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
    temp=adata_RNA.X[np.array(label)==celltype,:].mean(axis=0)
    TG=pd.DataFrame(temp.T,index=adata_RNA.var['gene_ids'].values,columns=['values'])
    del temp
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    TFbinding=load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN)
    if method=='LINGER':
        other_RE=list(set(N_overlap)-set(mat.index))
        if len(other_RE)>0:
            B_arr = pd.DataFrame(np.zeros((len(other_RE), mat.shape[1])), columns=mat.columns, index=other_RE)
            mat = pd.concat([mat, B_arr])
        mat = mat.loc[N_overlap]    
    if method=='baseline':
        mat=load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u)
        mat.index=N_overlap
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    TFbinding.index=N_overlap
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    TF_cluster = TF.values#[:,np.array(label)==celltype].mean(axis=1)
    TF_cluster=TF_cluster/TF_cluster.mean()
    #TF_cluster = TF_cluster[None,:]
    RE_cluster = RE.values#[:,np.array(label)==celltype].mean(axis=1)
    RE_cluster=RE_cluster/RE_cluster.mean()
    #RE_cluster = RE_cluster[:,None]
    S = (np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1)).T + np.log(TF_cluster+0.1)
    S = np.exp(S.T)
    S.index=N_overlap
    S_all = S.groupby(S.index).max()
    return S_all



def cell_type_specific_TF_RE_binding_score_scNN(mat,TFbinding,RE,TG,TFoverlap):
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    TF_cluster = TF.values#[:,np.array(label)==celltype].mean(axis=1)
    TF_cluster=TF_cluster/TF_cluster.mean()
    #TF_cluster = TF_cluster[None,:]
    RE_cluster = RE.values#[:,np.array(label)==celltype].mean(axis=1)
    RE_cluster=RE_cluster/RE_cluster.mean()
    #RE_cluster = RE_cluster[:,None]
    S = (np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1)).T + np.log(TF_cluster+0.1)
    S = np.exp(S.T)
    S.index=mat.index
    return S

def cell_type_specific_TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method):
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    if (celltype == 'all')&(method!='scNN'):
        for label0 in labelset:
            print('Generate cell type specitic TF binding potential for cell type '+ str(label0)+'...')
            result=pd.DataFrame()
            from tqdm import tqdm
            for i in tqdm(range(22)):
                chrN='chr'+str(i+1)
                mat=pd.read_csv(outdir+chrN+'_cell_population_TF_RE_binding.txt',sep='\t',index_col=0,header=0)
                out=cell_type_specific_TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,label0,outdir,method,mat)
        #result=pd.concat([result,out],axis=1).fillna(0)
                result = pd.concat([result, out], join='outer', axis=0)
            chrN='chrX'
            mat=pd.read_csv(outdir+chrN+'_cell_population_TF_RE_binding.txt',sep='\t',index_col=0,header=0)
            out=cell_type_specific_TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,label0,outdir,method,mat)
            result = pd.concat([result, out], join='outer', axis=0).fillna(0)
            result.to_csv(outdir+'cell_type_specific_TF_RE_binding_'+str(label0)+'.txt', sep='\t')
    elif method!='scNN':
        result=pd.DataFrame()
        from tqdm import tqdm
        chrom=['chr'+str(i+1) for i in range(22)]
        chrom.append('chrX')
        for i in tqdm(range(23)):
            chrN=chrom[i]
            mat=pd.read_csv(outdir+chrN+'_cell_population_TF_RE_binding.txt',sep='\t',index_col=0,header=0)
            out=cell_type_specific_TF_RE_binding_chr(adata_RNA,adata_ATAC,GRNdir,chrN,genome,celltype,outdir,method,mat)
        #result=pd.concat([result,out],axis=1).fillna(0)
            result = pd.concat([result, out], join='outer', axis=0)
        result.to_csv(outdir+'cell_type_specific_TF_RE_binding_'+str(celltype)+'.txt', sep='\t')
    elif (celltype == 'all')&(method=='scNN'):
        A=pd.read_csv(outdir+'cell_population_TF_RE_binding.txt',sep='\t',header=0,index_col=0)
        mat,REs,TFs=list2mat(A,'RE','TF','score')
        mat=pd.DataFrame(mat,index=REs,columns=TFs)
        TFs = mat.columns
        TFbinding=load_TFbinding_scNN(GRNdir,outdir,genome)
        TG=pd.DataFrame([],index=adata_RNA.var['gene_ids'].values)
        TFoverlap = list(set(TFs) & set(TG.index))
        TFoverlap=list(set(TFoverlap) & set(TFbinding.columns))
        mat = mat[TFoverlap]
        TFbinding = TFbinding[TFoverlap]
        REoverlap=list(set(TFbinding.index)&set(mat.index))
        TFbinding=TFbinding.loc[REoverlap]
        TFbinding1=np.zeros((mat.shape[0],len(TFoverlap)))
        REidx=pd.DataFrame(range(mat.shape[0]),index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values,:]=TFbinding.values
        TFbinding1 = pd.DataFrame(TFbinding1,index=mat.index,columns=TFoverlap)
        TFbinding=TFbinding1.copy()
        for label0 in labelset:
            print('Generate cell type specitic TF binding potential for cell type '+ str(label0)+'...')
            from tqdm import tqdm
            temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
            RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
            temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
            TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
            RE=RE.loc[REs]
            result=cell_type_specific_TF_RE_binding_score_scNN(mat,TFbinding,RE,TG,TFoverlap)
            result.to_csv(outdir+'cell_type_specific_TF_RE_binding_'+str(label0)+'.txt', sep='\t')
    else:
        label0=celltype
        A=pd.read_csv(outdir+'cell_population_TF_RE_binding.txt',sep='\t',header=0,index_col=0)
        mat,REs,TFs=list2mat(A,'RE','TF','score')
        mat=pd.DataFrame(mat,index=REs,columns=TFs)
        TFs = mat.columns
        TFbinding=load_TFbinding_scNN(GRNdir,outdir,genome)
        TG=pd.DataFrame([],index=adata_RNA.var['gene_ids'].values)
        TFoverlap = list(set(TFs) & set(TG.index))
        TFoverlap=list(set(TFoverlap) & set(TFbinding.columns))
        mat = mat[TFoverlap]
        TFbinding = TFbinding[TFoverlap]
        REoverlap=list(set(TFbinding.index)&set(RE.index))
        TFbinding=TFbinding.loc[REoverlap]
        TFbinding1=np.zeros((mat.shape[0],len(TFoverlap)))
        REidx=pd.DataFrame(range(mat.shape[0]),index=mat.index)
        TFbinding1[REidx.loc[TFbinding.index][0].values,:]=TFbinding.values
        TFbinding1 = pd.DataFrame(TFbinding1,index=mat.index,columns=TFoverlap)
        print('Generate cell type specitic TF binding potential for cell type '+ str(label0)+'...')
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        RE=RE.loc[REs]
        result=cell_type_specific_TF_RE_binding_score_scNN(mat,TFbinding,RE,TG,TFoverlap)
        result.to_csv(outdir+'cell_type_specific_TF_RE_binding_'+str(label0)+'.txt', sep='\t')
    

def load_shap(chr,outdir):
    import torch
    import pandas as pd
    import numpy as np
    import csv
    #print('loading shapley value '+chr+' ...')
    shap_all=torch.load(outdir+"shap_"+chr+".pt")
    import pandas as pd
    idx_file=outdir+'index.txt'
    TFName=outdir+'TFName.txt'
    #TFE=Input_dir+'TFexp.txt'
    REName = 'data/Peaks.txt'
# Open the file in read mode
    with open(REName, "r") as file:
    # Create a CSV reader
        reader = csv.reader(file, delimiter='\t')
    # Read the first column and store it in a list
        first_column = [row[0] for row in reader]
    REName=np.array(first_column)
    idx=pd.read_csv(idx_file,sep='\t',header=None)
    idx.fillna('', inplace=True)
    TFName=pd.read_csv(TFName,sep='\t',header=None)
    import numpy as np
    idx.columns=['gene','REid','TF_id','REid_b']
    TFName.columns=['Name']
    TFName=TFName['Name'].values
    #TFE=pd.read_csv(TFE,header=None,sep='\t')
    from scipy.stats import zscore
    TFindex=idx['TF_id'].values
    REindex=idx['REid'].values
    geneName=idx['gene'].values
    data_merge=pd.read_csv(outdir+'data_merge.txt',sep='\t',header=0,index_col=0)
    data_merge_temp=data_merge[data_merge['chr']==chr]
    return data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName
def cis_shap(chr,outdir):
    RE_2=[]
    TG_2=[]
    score_2=[]
    data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName=load_shap(chr,outdir)
    from tqdm import tqdm
    for j in tqdm(range(data_merge_temp.shape[0])):
        ii=data_merge_temp.index[j]
        if ii in shap_all.keys():
            AA0=shap_all[ii]
            REidxtemp=REindex[ii]
            REidxtemp=str(REidxtemp).split('_')
    #AA0[:,0:len(TFidxtemp)]=np.multiply(AA0[:,0:len(TFidxtemp)],TFE.values[np.array(TFidxtemp),:].T)
            temps=np.abs(AA0).mean(axis=0)
    #zscored_arr = zscore(temps)
            zscored_arr = np.nan_to_num(temps, nan=0.0)
            if (REidxtemp[0]=='') :
                REidxtemp=[]
            else:
                REidxtemp=[int(REidxtemp[i]) for i in range(len(REidxtemp))] 
            if len(REidxtemp)>0:
                REName_temp=REName[np.array(REidxtemp)]
                for k in range(len(REidxtemp)):
                    TG_2.append(geneName[ii])
                    RE_2.append(REName_temp[k])
                    score_2.append(zscored_arr[k+len(zscored_arr)-len(REidxtemp)])
    RE_TG=pd.DataFrame(TG_2)
    RE_TG.columns=['TG']
    RE_TG['RE']=RE_2
    RE_TG['score']=score_2
    RE_TG=RE_TG.groupby(['RE', 'TG'])['score'].max().reset_index()
    return RE_TG
def trans_shap(chr,outdir):
    TG_1=[]
    TF_1=[]
    score_1=[]
    data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName=load_shap(chr,outdir)
    from tqdm import tqdm
    for j in tqdm(range(data_merge_temp.shape[0])):
        ii=data_merge_temp.index[j]
        if ii in shap_all.keys():
            AA0=shap_all[ii]
            TFidxtemp=TFindex[ii]
            TFidxtemp=TFidxtemp.split('_')
            TFidxtemp=[int(TFidxtemp[i]) for i in range(len(TFidxtemp))]
            TFName_temp=TFName[np.array(TFidxtemp)]
    #AA0[:,0:len(TFidxtemp)]=np.multiply(AA0[:,0:len(TFidxtemp)],TFE.values[np.array(TFidxtemp),:].T)
            temps=np.abs(AA0).mean(axis=0)
    #zscored_arr = zscore(temps)
            zscored_arr = np.nan_to_num(temps, nan=0.0)
        for k in range(len(TFidxtemp)):
            TG_1.append(geneName[ii])
            TF_1.append(TFName_temp[k])
            score_1.append(zscored_arr[k])
    TF_TG=pd.DataFrame(TG_1)
    TF_TG.columns=['TG']
    TF_TG['TF']=TF_1
    TF_TG['score']=score_1
    mat,TGs,TFs=list2mat(TF_TG,'TG','TF','score')
    mat=pd.DataFrame(mat,index=TGs,columns=TFs)
    mat.fillna(0, inplace=True)
    return mat
      
def load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap):   
    #print('load prior RE-TG ...')
    from scipy.sparse import coo_matrix
    primary_s=pd.read_csv(GRNdir+'Primary_RE_TG_'+chrN+'.txt',sep='\t')
    primary_s["RE"] = primary_s["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    primary_s = primary_s[primary_s["RE"].isin(O_overlap_u)]
    TGset=primary_s["TG"].unique()
    REset=O_overlap_u
    # Create a dictionary mapping column names and row names to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
# Map the column names and row names to integer indices in the DataFrame
    primary_s.loc[:,"col_index"] = primary_s["TG"].map(col_dict)
    primary_s.loc[:,"row_index"] = primary_s["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = primary_s["col_index"].tolist()
    row_indices = primary_s["row_index"].tolist()
    values = primary_s["score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
    array = sparse_S.toarray()
    O_overlap_u_df=pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38=pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    array2=np.zeros([len(O_overlap),array.shape[1]])
    index=O_overlap_u_df.loc[O_overlap][0].values
    array2=array[index,:]
    array=pd.DataFrame(array2,index=O_overlap,columns=TGset)
    return array,TGset
def load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGoverlap):
    #print('load RE-TG distance for '+chrN+'...')
    from scipy.sparse import coo_matrix
    Dis=pd.read_csv(GRNdir+'RE_TG_distance_'+chrN+'.txt',sep='\t',header=None)
    Dis.columns=['RE','TG','dis']
    Dis["RE"] = Dis["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    Dis = Dis[Dis["RE"].isin(O_overlap_hg19_u)]
    Dis = Dis[Dis['TG'].isin(TGoverlap)]
    col_dict = {col: i for i, col in enumerate(TGoverlap)}
    row_dict = {row: i for i, row in enumerate(O_overlap_hg19_u)}
# Map the column names and row names to integer indices in the DataFrame
    Dis.loc[:,"col_index"] = Dis["TG"].map(col_dict)
    Dis.loc[:,"row_index"] = Dis["RE"].map(row_dict)
    col_indices = Dis["col_index"].tolist()
    row_indices = Dis["row_index"].tolist()
    values = Dis["dis"].tolist()
# Create the sparse matrix using coo_matrix
    sparse_dis = coo_matrix((values, (row_indices, col_indices)),shape=(len(O_overlap_u), len(TGoverlap)))
    sparse_dis.colnames = TGoverlap
    sparse_dis.rownames = O_overlap_u
    sparse_dis = sparse_dis.tocsc()
    A=sparse_dis.multiply(1 / 25000)
    A.data +=0.5
    A.data = np.exp(-A.data)
    sparse_dis=A
    array = sparse_dis.toarray()
    O_overlap_u_df=pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38=pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    array2=np.zeros([len(O_overlap),array.shape[1]])
    index=O_overlap_u_df.loc[O_overlap][0].values
    array2=array[index,:]
    array=pd.DataFrame(array2,index=O_overlap,columns=TGoverlap)
    return array


def load_RE_TG_scNN(outdir):   
    #print('load prior RE-TG ...')
    from scipy.sparse import coo_matrix
    import pandas as pd
    import numpy as np
    dis=pd.read_csv('data/RE_gene_distance.txt',sep='\t',header=0)
    dis['distance']=np.exp(-(0.5+dis['distance']/25000))
    REs=dis['RE'].unique()
    TGs=dis['gene'].unique()
    cis=pd.read_csv(outdir+'cell_population_cis_regulatory.txt',sep='\t',header=None)
    cis.columns=['RE','TG','score']
    REs2=cis['RE'].unique()
    TGs2=cis['TG'].unique()
    REoverlap=list(set(REs2)&set(REs))
    TGoverlap=list(set(TGs)&set(TGs2))
    cisGRN,REs2,TGs2=list2mat_s(cis,REoverlap,TGoverlap,'RE','TG','score')
    dis=dis[dis['RE'].isin(REoverlap)]
    dis=dis[dis['gene'].isin(TGoverlap)]
    distance,REs,TGs=list2mat_s(dis,REoverlap,TGoverlap,'RE','gene','distance')
    return distance,cisGRN,REoverlap,TGoverlap

def cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,outdir):  
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(GRNdir,genome,chrN,outdir)
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    RE=pd.DataFrame(adata_ATAC.X.toarray().T,index=adata_ATAC.var['gene_ids'].values,columns=adata_ATAC.obs['barcode'].values)
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    RE=RE.mean(axis=1)
    RE=RE/RE.mean()+0.1
    ## select the genes
    TG=pd.DataFrame(adata_RNA.X.toarray().T,index=adata_RNA.var['gene_ids'].values,columns=adata_RNA.obs['barcode'].values)
    TGoverlap=list(set(TGset)&set(TG.index))
    #target_col_indices = [col_dict[col] for col in TGoverlap]
    sparse_S = sparse_S[TGoverlap]
    TG=TG.loc[TGoverlap]
    TG=TG.mean(axis=1)
    TG=TG/TG.mean()+0.1
    sparse_dis=load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGoverlap)
    sparse_S+=0.1
    #Score=csc_matrix(RE).T.multiply(sparse_S.values).multiply(sparse_dis.values).multiply(csc_matrix(TG)).toarray()
    Score=np.multiply(sparse_S.values,sparse_dis.values)
    Score=pd.DataFrame(Score,index=N_overlap,columns=TGoverlap)
    Score=Score.groupby(Score.index).max()
    data = Score.values[Score.values!=0] 
    rows, cols = np.nonzero(Score.values) 
    coo = coo_matrix((data,(rows,cols)),shape=Score.shape)
    combined = np.zeros([len(data),3], dtype=object) 
    combined[:,0]=Score.index[coo.row]
    combined[:,1]=np.array(TGoverlap)[coo.col]
    combined[:,2]=coo.data
    combined=pd.DataFrame(combined)
    return combined

def cis_shap_scNN(chrtemp,outdir,RE_TGlink1,REName,TFName):
    import ast
    REName=pd.DataFrame(range(len(REName)),index=REName)
    RE_2=[]
    TG_2=[]
    score_2=[]
    shap_all=torch.load(outdir+chrtemp+"_shap"+".pt")
    N=RE_TGlink1.shape[0]
    for ii in tqdm(range(N)):
        AA0=shap_all[ii]
        RE_TGlink_temp=RE_TGlink1.values[ii,:]
        actual_list = ast.literal_eval(RE_TGlink_temp[1])
        REidxtemp=REName.loc[actual_list].index
        TFidxtemp=np.array(range(len(TFName)))
        TFidxtemp=TFidxtemp[TFName!=RE_TGlink_temp[0]]
        if len(REidxtemp)>0:
            temps=np.abs(AA0).mean(axis=0)
            zscored_arr = np.nan_to_num(temps, nan=0.0)
            for k in range(len(REidxtemp)):
                TG_2.append(RE_TGlink_temp[0])
                RE_2.append(REidxtemp[k])
                score_2.append(zscored_arr[k+len(zscored_arr)-len(REidxtemp)])
    RE_TG=pd.DataFrame(TG_2)
    RE_TG.columns=['TG']
    RE_TG['RE']=RE_2
    RE_TG['score']=score_2
    RE_TG=RE_TG.groupby(['RE', 'TG'])['score'].max().reset_index()
    return RE_TG


def cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir): 
    from tqdm import tqdm
    chrom=['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    if method=='baseline':
        result=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            temp=cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,outdir)
            temp.columns=['RE','TG','Score']
            result=pd.concat([result,temp],axis=0,join='outer')
    if method=='LINGER':
        result=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            temp=cis_shap(chrN,outdir)
            result=pd.concat([result,temp],axis=0,join='outer')
    if method=='scNN':
        Exp,Opn,Target,RE_TGlink=load_data_scNN(GRNdir,genome)
        RE_TGlink=pd.read_csv(outdir+'RE_TGlink.txt',sep='\t',header=0)
        RE_TGlink.columns=[0,1,'chr']
        #chrall=[RE_TGlink[0][i][0].split(':')[0] for i in range(RE_TGlink.shape[0])]
        chrlist=RE_TGlink['chr'].unique()
        REName=Opn.index
        geneName=Target.index
        TFName=Exp.index
        result=pd.DataFrame([])
        for i in tqdm(range(len(chrlist))):
            chrN=chrlist[i]
            RE_TGlink1=RE_TGlink[RE_TGlink['chr']==chrN]
            temp=cis_shap_scNN(chrN,outdir,RE_TGlink1,REName,TFName)
            result=pd.concat([result,temp],axis=0,join='outer')
    result.to_csv(outdir+'cell_population_cis_regulatory.txt',sep='\t',header=None,index=None)


def cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,celltype,outdir): 
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(GRNdir,genome,chrN,outdir)
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    temp=adata_ATAC.X[np.array(label)==celltype,:].mean(axis=0).T
    RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
    temp=adata_RNA.X[np.array(label)==celltype,:].mean(axis=0).T
    TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
    del temp
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    ## select the genes
    TGoverlap=list(set(TGset)&set(TG.index))
    #target_col_indices = [col_dict[col] for col in TGoverlap]
    sparse_S = sparse_S[TGoverlap]
    TG=TG.loc[TGoverlap]
    sparse_dis=load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGoverlap)
    sparse_S+=0.1
## cell annotation
    TG_temp=TG.values#[:,np.array(label)==celltype].mean(axis=1)
    TG_temp=TG_temp/TG_temp.mean()+0.1
    RE_temp=RE.values#[:,np.array(label)==celltype].mean(axis=1)
    RE_temp=RE_temp/RE_temp.mean()+0.1
    Score=csc_matrix(RE_temp).multiply(sparse_S.values).multiply(sparse_dis.values).multiply(csc_matrix(TG_temp.T)).toarray()
    Score=pd.DataFrame(Score,index=N_overlap,columns=TGoverlap)
    Score=Score.groupby(Score.index).max()
    data = Score.values[Score.values!=0] 
    rows, cols = np.nonzero(Score.values) 
    coo = coo_matrix((data,(rows,cols)),shape=Score.shape)
    combined = np.zeros([len(data),3], dtype=object) 
    combined[:,0]=Score.index[coo.row]
    combined[:,1]=np.array(TGoverlap)[coo.col]
    combined[:,2]=coo.data
    resultall=pd.DataFrame(combined)
    return resultall  
def cell_type_specific_cis_reg_scNN(distance,cisGRN,RE,TG,REs,TGs):
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    from scipy.sparse import coo_matrix
    RE=RE.loc[REs]
    ## select the genes
    #target_col_indices = [col_dict[col] for col in TGoverlap]
    TG=TG.loc[TGs]
    ## cell annotation
    TG_temp=TG.values#[:,np.array(label)==celltype].mean(axis=1)
    TG_temp=TG_temp/TG_temp.mean()+0.1
    RE_temp=RE.values#[:,np.array(label)==celltype].mean(axis=1)
    RE_temp=RE_temp/RE_temp.mean()+0.1
    Score=(cisGRN.multiply(csr_matrix(RE_temp))).multiply(distance).multiply(csr_matrix(TG_temp.T))
    row_indices, col_indices = Score.nonzero()
    row_indices=np.array(REs)[row_indices]
    col_indices = np.array(TGs)[col_indices]
    values = Score.data
    combined = np.zeros([len(row_indices),3], dtype=object) 
    combined[:,0]=row_indices
    combined[:,1]=col_indices
    combined[:,2]=values
    resultall=pd.DataFrame(combined)
    return resultall 

def cell_type_specific_cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method): 
    import pandas as pd
    import numpy as np
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    chrom=['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    from tqdm import tqdm
    if (celltype=='all')&(method!='scNN'):
        for label0 in labelset:
            label0=str(label0)
            result=pd.DataFrame([])
            for i in tqdm(range(23)):
                chrN=chrom[i]
                temp=cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,label0,outdir)
                result=pd.concat([result,temp],axis=0,join='outer')
            chrN='chrX'
            temp=cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,label0,outdir)
            result=pd.concat([result,temp],axis=0,join='outer')
            result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+str(label0)+'.txt',sep='\t',header=None,index=None)
    elif (method!='scNN'):
            result=pd.DataFrame([])
            for i in tqdm(range(23)):
                chrN=chrom[i]
                temp=cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,celltype,outdir)
                result=pd.concat([result,temp],axis=0,join='outer')
            chrN='chrX'
            temp=cell_type_specific_cis_reg_chr(GRNdir,adata_RNA,adata_ATAC,genome,chrN,celltype,outdir)
            result=pd.concat([result,temp],axis=0,join='outer')
            result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+celltype+'.txt',sep='\t',header=None,index=None)
    elif (celltype=='all')&(method=='scNN'):
        distance,cisGRN,REs,TGs=load_RE_TG_scNN(outdir)
        for label0 in labelset:
            label0=str(label0)
            temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
            RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
            temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
            TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
            del temp
            
            result=cell_type_specific_cis_reg_scNN(distance,cisGRN,RE,TG,REs,TGs)
            result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+label0+'.txt',sep='\t',header=None,index=None)
    else: 
        label0=celltype
        label0=str(label0)
        temp=adata_ATAC.X[np.array(label)==label0,:].mean(axis=0).T
        RE=pd.DataFrame(temp,index=adata_ATAC.var['gene_ids'].values,columns=['values'])
        temp=adata_RNA.X[np.array(label)==label0,:].mean(axis=0).T
        TG=pd.DataFrame(temp,index=adata_RNA.var['gene_ids'].values,columns=['values'])
        del temp
        result=cell_type_specific_cis_reg_scNN(distance,cisGRN,RE,TG,REs,TGs)
        result.to_csv(outdir+'cell_type_specific_cis_regulatory_'+label0+'.txt',sep='\t',header=None,index=None)
            
def trans_shap_scNN(chrtemp,outdir,RE_TGlink1,REName,TFName):
    import ast
    TG_1=[]
    TF_1=[]
    score_1=[]
    REName=pd.DataFrame(range(len(REName)),index=REName)
    shap_all=torch.load(outdir+chrtemp+"_shap"+".pt")
    N=RE_TGlink1.shape[0]
    from tqdm import tqdm
    for ii in range(N):
        AA0=shap_all[ii]
        RE_TGlink_temp=RE_TGlink1.values[ii,:]
        actual_list = ast.literal_eval(RE_TGlink_temp[1])
        REidxtemp=REName.loc[actual_list].index
        TFidxtemp=np.array(range(len(TFName)))
        TFidxtemp=TFidxtemp[TFName!=RE_TGlink_temp[0]]
        temps=np.abs(AA0).mean(axis=0)
    #zscored_arr = zscore(temps)
        zscored_arr = np.nan_to_num(temps, nan=0.0)
        for k in range(len(TFidxtemp)):
            TG_1.append(RE_TGlink_temp[0])
            TF_1.append(TFName[TFidxtemp[k]])
            score_1.append(zscored_arr[k])
    TF_TG=pd.DataFrame(TG_1)
    TF_TG.columns=['TG']
    TF_TG['TF']=TF_1
    TF_TG['score']=score_1
    mat,TGs,TFs=list2mat(TF_TG,'TG','TF','score')
    mat=pd.DataFrame(mat,index=TGs,columns=TFs)
    mat.fillna(0, inplace=True)
    return mat


def load_cis(Binding,celltype,outdir):
    from scipy.sparse import coo_matrix
    import pandas as pd
    import numpy as np
    if celltype=='':
        cis=pd.read_csv(outdir+'cell_population_cis_regulatory.txt',sep='\t',header=None)
    else:
        cis=pd.read_csv(outdir+'cell_type_specific_cis_regulatory_'+celltype+'.txt',sep='\t',header=None)
    cis.columns=['RE','TG','Score']
    TGset=cis['TG'].unique()
    REset=Binding.index
    TFset=Binding.columns
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
    cis=cis[cis["RE"].isin(REset)]
# Map the column names and row names to integer indices in the DataFrame
    cis["col_index"] = cis["TG"].map(col_dict)
    cis["row_index"] = cis["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = cis["col_index"].tolist()
    row_indices = cis["row_index"].tolist()
    values = cis["Score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)),shape=(len(REset), len(TGset)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
    cis=sparse_S.toarray()
    cis=pd.DataFrame(cis,index=REset,columns=TGset)
    return cis

def load_TF_TG( GRNdir, TFset,TGset):
    TF_TG_all=np.zeros([len(TGset),len(TFset)])
    a=list(range(1,23))
    a.append('X')
    for i in a:
        chrN='chr'+str(i)     
        TF_TG = pd.read_csv(GRNdir+'Primary_TF_TG_'+chrN+'.txt',sep='\t')
        TF_TG = TF_TG[TF_TG['TF'].isin(TFset)]
        TF_TG = TF_TG[TF_TG['TG'].isin(TGset)]
        col_dict = {col: i for i, col in enumerate(TFset)}
        row_dict = {row: i for i, row in enumerate(TGset)}
        TF_TG["col_index"] = TF_TG["TF"].map(col_dict)
        TF_TG["row_index"] = TF_TG["TG"].map(row_dict)
        col_indices = TF_TG["col_index"].tolist()
        row_indices = TF_TG["row_index"].tolist()
        values = TF_TG["score"].tolist()
        sparse_S = coo_matrix((values, (row_indices, col_indices)),shape=( len(TGset),len(TFset)))
        idx=list(set(row_indices))
        TGset1=TGset[idx]
        TF_TG=sparse_S.toarray()
        TF_TG=pd.DataFrame(TF_TG,index=TGset,columns=TFset)
        TF_TG=TF_TG.loc[TGset1]
        TF_TG_all[idx,:]=TF_TG.values
    TF_TG_all=pd.DataFrame(TF_TG_all,index=TGset,columns=TFset)
    return TF_TG_all

def trans_reg(GRNdir,method,outdir,genome):
    import ast
    import pandas as pd
    from scipy.sparse import coo_matrix
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    print('Generate trans-regulatory netowrk ...')
    if method=='baseline':
        Binding=pd.read_csv(outdir+'cell_population_TF_RE_binding.txt',sep='\t',index_col=0)
        cis=load_cis(Binding,'',outdir)
        TFset=Binding.columns
        TGset=cis.columns
        TF_TG=load_TF_TG(GRNdir, TFset,TGset)
        S=np.matmul(Binding.values.T, cis.values).T*(TF_TG.values.T).T
        S=pd.DataFrame(S, index=TGset,columns=TFset)
    elif method=='LINGER':
        chrom=['chr'+str(i+1) for i in range(22)]
        chrom.append('chrX')
        S=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            temp=trans_shap(chrN,outdir)
            S=pd.concat([S,temp],axis=0,join='outer')
    elif method=='scNN':
        Exp,Opn,Target,RE_TGlink=load_data_scNN(GRNdir,genome)
        RE_TGlink=pd.read_csv(outdir+'RE_TGlink.txt',sep='\t',header=0)
        RE_TGlink.columns=[0,1,'chr']
        #chrall=[RE_TGlink[0][i][0].split(':')[0] for i in range(RE_TGlink.shape[0])]
        chrlist=RE_TGlink['chr'].unique()
        REName=Opn.index
        geneName=Target.index
        TFName=Exp.index
        result=pd.DataFrame([])
        S=pd.DataFrame([])
        for i in tqdm(range(len(chrlist))):
            chrN=chrlist[i]
            RE_TGlink1=RE_TGlink[RE_TGlink['chr']==chrN]
            temp=trans_shap_scNN(chrN,outdir,RE_TGlink1,REName,TFName)
            S=pd.concat([S,temp],axis=0,join='outer')
    print('Save trans-regulatory netowrk ...')
    S.to_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t')

def cell_type_specific_trans_reg(GRNdir,adata_RNA,celltype,outdir):
    import pandas as pd
    import numpy as np
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    label=adata_RNA.obs['label'].values.tolist()
    labelset=list(set(label))
    if celltype=='all':
        for label0 in labelset:
            Binding=pd.read_csv(outdir+'cell_type_specific_TF_RE_binding_'+str(label0)+'.txt',sep='\t',index_col=0)
            label0=str(label0)
            cis=load_cis(Binding,label0,outdir)
            TFset=Binding.columns
            TGset=cis.columns
            #TF_TG=load_TF_TG(GRNdir, TFset,TGset)
            S=np.matmul(Binding.values.T, cis.values).T#*(TF_TG.values.T).T
            S=pd.DataFrame(S, index=TGset,columns=TFset)
            S.to_csv(outdir+'cell_type_specific_trans_regulatory_'+str(label0)+'.txt',sep='\t')
    else:
        Binding=pd.read_csv(outdir+'cell_type_specific_TF_RE_binding_'+celltype+'.txt',sep='\t',index_col=0)
        cis=load_cis(Binding,celltype,outdir)
        TFset=Binding.columns
        TGset=cis.columns
        #TF_TG=load_TF_TG(GRNdir, TFset,TGset)
        S=np.matmul(Binding.values.T, cis.values).T#*(TF_TG.values.T).T
        S=pd.DataFrame(S, index=TGset,columns=TFset)
        S.to_csv(outdir+'cell_type_specific_trans_regulatory_'+celltype+'.txt',sep='\t')

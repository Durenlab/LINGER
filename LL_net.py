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
from sklearn.metrics.pairwise import cosine_similarity
hidden_size  = 64
hidden_size2 = 16
output_size = 1
seed_value = 42
torch.manual_seed(seed_value)
class Net(nn.Module):
    def __init__(self,input_size):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.fc2 = nn.Linear(64, 16)
        self.fc3 = nn.Linear(16, output_size)
    def forward(self, x):
        #x = torch.sigmoid(self.fc1(x))
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x
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
        corr_matrix = cosine_similarity(temps.detach().numpy().T)
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

def load_region(Input_dir,GRNdir,genome,chrN):  
    O_overlap=merge_columns_in_bed_file(Input_dir+'Region_overlap_'+chrN+'.bed',1)
    N_overlap=merge_columns_in_bed_file(Input_dir+'Region_overlap_'+chrN+'.bed',4)
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
def TF_RE_LINGER_chr(Input_dir,ATAC_file,chr):
    REName = Input_dir+ATAC_file
# Open the file in read mode
    with open(REName, "r") as file:
    # Create a CSV reader
        reader = csv.reader(file, delimiter='\t')
    # Read the first column and store it in a list
        first_column = [row[0] for row in reader]
    REName=np.array(first_column)
    idx_file=Input_dir+'index.txt'
    from scipy.stats import zscore
    data0=pd.read_csv(Input_dir+'result_'+chr+'.txt',sep='\t')
    data0.columns=['gene','x','y']
    idx_file=Input_dir+'index.txt'
    idx=pd.read_csv(idx_file,sep='\t',header=None)
    idx.columns=['gene','REid','TF_id','REid_b']
    idx.fillna('', inplace=True)
    TFName=Input_dir+'TFName.txt'
    TFName=pd.read_csv(TFName,sep='\t',header=None)
    TFName.columns=['Name']
    TFName=TFName['Name'].values
    TFindex=idx['TF_id'].values
    REindex=idx['REid'].values
    geneName=idx['gene'].values
    net_all=torch.load(Input_dir+"net_"+chr+".pt")
    data_merge=pd.read_csv(Input_dir+'data_merge.txt',sep='\t',header=0,index_col=0)
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

def TF_RE_binding_chr(RNA_file,ATAC_file,Input_dir,GRNdir,chrN,genome):
    ## the regions
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(Input_dir,GRNdir,genome,chrN)
    import numpy as np
    import pandas as pd
## read the count file.
    RE=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    TFbinding=load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN)
    mat=load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u)
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    TF_cluster = TF.values.mean(axis=1)
    TF_cluster = TF_cluster[None,:]
    RE_cluster = RE.values.mean(axis=1)
    RE_cluster = RE_cluster[:,None]
    S = np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1) + np.log(TF_cluster+0.1)
    S = np.exp(S)
    S.index=N_overlap
    mean_S = S.groupby(S.index).max()
    return mean_S
def TF_RE_binding(Input_dir,GRNdir,RNA_file,ATAC_file,genome,method):
    from tqdm import tqdm
    print('Generating cellular population TF binding strength ...')
    chrom = ['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    if method=='baseline':
        result=pd.DataFrame()
        for i in tqdm(range(23)):
            chrN=chrom[i]
            out=TF_RE_binding_chr(RNA_file,ATAC_file,Input_dir,GRNdir,chrN,genome)
        #result=pd.concat([result,out],axis=1).fillna(0)
            result = pd.concat([result, out], join='outer', axis=0)
        print('Generating cellular population TF binding strength for chrX')
        chrN='chrX'
        out=TF_RE_binding_chr(RNA_file,ATAC_file,Input_dir,GRNdir,chrN,genome)
        result = pd.concat([result, out], join='outer', axis=0)
    if method=='LINGER':
        result=pd.DataFrame()
        for i in tqdm(range(23)):
            chr=chrom[i]
            mat=TF_RE_LINGER_chr(Input_dir,ATAC_file,chr)
            result = pd.concat([result, mat], join='outer', axis=0)
    result.to_csv(Input_dir+'cell_population_TF_RE_binding.txt',sep='\t')


def cell_type_specific_TF_RE_binding_chr(RNA_file,ATAC_file,labels,Input_dir,GRNdir,chrN,genome,celltype):
    ## the regions
        ## the regions
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(Input_dir,GRNdir,genome,chrN)
    import numpy as np
    import pandas as pd
## read the count file.
    RE=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    TFbinding=load_TFbinding(GRNdir,O_overlap,O_overlap_u,O_overlap_hg19_u,chrN)
    mat=load_TF_RE(GRNdir,chrN,O_overlap,O_overlap_u,O_overlap_hg19_u)
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    with open(Input_dir+labels, 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
    TF_cluster = TF.values[:,np.array(label)==celltype].mean(axis=1)
    TF_cluster=TF_cluster/TF_cluster.mean()
    TF_cluster = TF_cluster[None,:]
    RE_cluster = RE.values[:,np.array(label)==celltype].mean(axis=1)
    RE_cluster=RE_cluster/RE_cluster.mean()
    RE_cluster = RE_cluster[:,None]
    S = np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1) + np.log(TF_cluster+0.1)
    S = np.exp(S)
    S.index=N_overlap
    S_all = S.groupby(S.index).max()
    return S_all

def cell_type_specific_TF_RE_binding(Input_dir,GRNdir,RNA_file,ATAC_file,labels,genome,celltype):
    with open(Input_dir+labels, 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
    if celltype == 'all':
        for label0 in range(len(labelset)):
            print('Generate cell type specitic TF binding potential for cell type '+ label0+'...')
            result=pd.DataFrame()
            from tqdm import tqdm
            for i in tqdm(range(22)):
                chrN='chr'+str(i+1)
                out=cell_type_specific_TF_RE_binding_chr(RNA_file,ATAC_file,labels,Input_dir,GRNdir,chrN,genome,label0)
        #result=pd.concat([result,out],axis=1).fillna(0)
                result = pd.concat([result, out], join='outer', axis=0)
            chrN='chrX'
            out=cell_type_specific_TF_RE_binding_chr(RNA_file,ATAC_file,labels,Input_dir,GRNdir,chrN,genome,label0)
            result = pd.concat([result, out], join='outer', axis=0).fillna(0)
            result.to_csv(Input_dir+'cell_type_specific_TF_RE_binding_'+label0+'.txt', sep='\t')
    else:
        result=pd.DataFrame()
        from tqdm import tqdm
        chrom=['chr'+str(i+1) for i in range(22)]
        chrom.append('chrX')
        for i in tqdm(range(23)):
            chrN=chrom[i]
            out=cell_type_specific_TF_RE_binding_chr(RNA_file,ATAC_file,labels,Input_dir,GRNdir,chrN,genome,celltype)
        #result=pd.concat([result,out],axis=1).fillna(0)
            result = pd.concat([result, out], join='outer', axis=0)
        result.to_csv(Input_dir+'cell_type_specific_TF_RE_binding_'+celltype+'.txt', sep='\t')
    

def load_shap(chr,Input_dir):
    import torch
    import pandas as pd
    import numpy as np
    import csv
    #print('loading shapley value '+chr+' ...')
    shap_all=torch.load(Input_dir+"shap_"+chr+".pt")
    import pandas as pd
    idx_file=Input_dir+'index.txt'
    TFName=Input_dir+'TFName.txt'
    #TFE=Input_dir+'TFexp.txt'
    REName = Input_dir+"ATAC.txt"
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
    data_merge=pd.read_csv(Input_dir+'data_merge.txt',sep='\t',header=0,index_col=0)
    data_merge_temp=data_merge[data_merge['chr']==chr]
    return data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName
def cis_shap(Input_dir,chr):
    RE_2=[]
    TG_2=[]
    score_2=[]
    data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName=load_shap(chr,Input_dir)
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
def trans_shap(Input_dir,chr):
    TG_1=[]
    TF_1=[]
    score_1=[]
    data_merge_temp,geneName,REindex,TFindex,shap_all,TFName,REName=load_shap(chr,Input_dir)
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
def cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN):  
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(Input_dir,GRNdir,genome,chrN)
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    ## read the count file.
    RE=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
## cell annotation
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    ## select the genes
    TGoverlap=list(set(TGset)&set(TG.index))
    #target_col_indices = [col_dict[col] for col in TGoverlap]
    sparse_S = sparse_S[TGoverlap]
    TG=TG.loc[TGoverlap]
    sparse_dis=load_RE_TG_distance(GRNdir,chrN,O_overlap_hg19_u,O_overlap_u,O_overlap,TGoverlap)
    TG=TG.mean(axis=1)
    TG=TG/TG.mean()+0.1
    RE=RE.mean(axis=1)
    RE=RE/RE.mean()+0.1
    sparse_S+=0.1
    Score=csc_matrix(RE).T.multiply(sparse_S.values).multiply(sparse_dis.values).multiply(csc_matrix(TG)).toarray()
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
def cis_reg(Input_dir,GRNdir,RNA_file,ATAC_file,genome,method): 
    from tqdm import tqdm
    chrom=['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    if method=='baseline':
        result=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            temp=cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN)
            temp.columns=['RE','TG','Score']
            result=pd.concat([result,temp],axis=0,join='outer')
    if method=='LINGER':
        result=pd.DataFrame([])
        for i in tqdm(range(23)):
            chrN=chrom[i]
            temp=cis_shap(Input_dir,chrN)
            result=pd.concat([result,temp],axis=0,join='outer')
    result.to_csv(Input_dir+'cell_population_cis_regulatory.txt',sep='\t',header=None,index=None)


def cell_type_specific_cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN,celltype): 
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    O_overlap, N_overlap,O_overlap_u,N_overlap_u,O_overlap_hg19_u=load_region(Input_dir,GRNdir,genome,chrN)
    sparse_S,TGset=load_RE_TG(GRNdir,chrN,O_overlap_u,O_overlap_hg19_u,O_overlap)
    ## read the count file.
    RE=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
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
    with open(Input_dir+'label.txt', 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
    TG_temp=TG.values[:,np.array(label)==celltype].mean(axis=1)
    TG_temp=TG_temp/TG_temp.mean()+0.1
    RE_temp=RE.values[:,np.array(label)==celltype].mean(axis=1)
    RE_temp=RE_temp/RE_temp.mean()+0.1
    Score=csc_matrix(RE_temp).T.multiply(sparse_S.values).multiply(sparse_dis.values).multiply(csc_matrix(TG_temp)).toarray()
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

def cell_type_specific_cis_reg(Input_dir,GRNdir,RNA_file,ATAC_file,label_file,genome,celltype): 
    import pandas as pd
    import numpy as np
    with open(Input_dir+label_file, 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
    chrom=['chr'+str(i+1) for i in range(22)]
    chrom.append('chrX')
    from tqdm import tqdm
    if celltype=='all':
        for label0 in labelset:
            result=pd.DataFrame([])
            for i in tqdm(range(23)):
                chrN=chrom[i]
                temp=cell_type_specific_cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN,label0)
                result=pd.concat([result,temp],axis=0,join='outer')
            chrN='chrX'
            temp=cell_type_specific_cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN,label0)
            result=pd.concat([result,temp],axis=0,join='outer')
            result.to_csv(Input_dir+'cell_type_specific_cis_regulatory_'+label0+'.txt',sep='\t',header=None,index=None)
    else:
            result=pd.DataFrame([])
            for i in tqdm(range(23)):
                chrN=chrom[i]
                temp=cell_type_specific_cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN,celltype)
                result=pd.concat([result,temp],axis=0,join='outer')
            chrN='chrX'
            temp=cell_type_specific_cis_reg_chr(Input_dir,GRNdir,RNA_file,ATAC_file,genome,chrN,celltype)
            result=pd.concat([result,temp],axis=0,join='outer')
            result.to_csv(Input_dir+'cell_type_specific_cis_regulatory_'+celltype+'.txt',sep='\t',header=None,index=None)



def load_cis(Input_dir,Binding,celltype):
    if celltype=='':
        cis=pd.read_csv(Input_dir+'cell_population_cis_regulatory.txt',sep='\t',header=None)
    else:
        cis=pd.read_csv(Input_dir+'cell_type_specific_cis_regulatory_'+celltype+'.txt',sep='\t',header=None)
    cis.columns=['RE','TG','Score']
    TGset=cis['TG'].unique()
    REset=Binding.index
    TFset=Binding.columns
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
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

def trans_reg(Input_dir,GRNdir,RNA_file,ATAC_file,method):
    import pandas as pd
    from scipy.sparse import coo_matrix
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    print('Generate trans-regulatory netowrk ...')
    if method=='baseline':
        Binding=pd.read_csv(Input_dir+'cell_population_TF_RE_binding.txt',sep='\t',index_col=0)
        cis=load_cis(Input_dir,Binding,'')
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
            temp=trans_shap(Input_dir,chrN)
            S=pd.concat([S,temp],axis=0,join='outer')
    print('Save trans-regulatory netowrk ...')
    S.to_csv(Input_dir+'cell_population_trans_regulatory.txt',sep='\t')

def cell_type_specific_trans_reg(Input_dir,GRNdir,RNA_file,labels,ATAC_file,celltype):
    import pandas as pd
    from scipy.sparse import coo_matrix
    import numpy as np
    import pandas as pd
    from scipy.sparse import csc_matrix
    from scipy.sparse import coo_matrix
    with open(Input_dir+'label.txt', 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
    if celltype=='all':
        for label0 in labelset:
            Binding=pd.read_csv(Input_dir+'cell_type_specific_TF_RE_binding_'+label0+'.txt',sep='\t',index_col=0)
            cis=load_cis(Input_dir,Binding,label0)
            TFset=Binding.columns
            TGset=cis.columns
            TF_TG=load_TF_TG(GRNdir, TFset,TGset)
            S=np.matmul(Binding.values.T, cis.values).T*(TF_TG.values.T).T
            S=pd.DataFrame(S, index=TGset,columns=TFset)
            S.to_csv(Input_dir+'cell_type_specific_trans_regulatory_'+label0+'.txt',sep='\t')
    else:
        Binding=pd.read_csv(Input_dir+'cell_type_specific_TF_RE_binding_'+celltype+'.txt',sep='\t',index_col=0)
        cis=load_cis(Input_dir,Binding,celltype)
        TFset=Binding.columns
        TGset=cis.columns
        TF_TG=load_TF_TG(GRNdir, TFset,TGset)
        S=np.matmul(Binding.values.T, cis.values).T*(TF_TG.values.T).T
        S=pd.DataFrame(S, index=TGset,columns=TFset)
        S.to_csv(Input_dir+'cell_type_specific_trans_regulatory_'+celltype+'.txt',sep='\t')

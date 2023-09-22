import os
import numpy as np
import pandas as pd
from immupute_dis import immupute_dis
import pseudo_bulk
import subprocess
from tqdm import tqdm
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



def gene_expression(Input_dir,GRNdir,TFName,Match2,TG_psedubulk):
	gene = pd.read_csv(GRNdir+'bulk_gene_all.txt')
	gene.columns=['gene']
	#gene=gene['gene']
	d1 = np.isin(TG_psedubulk.index, gene['gene'].values)
	List = TG_psedubulk.index[d1]
	A = np.log2(1 + TG_psedubulk.values[d1, :])
	#Write Exp.txt and Symbol.txt
	pd.DataFrame(A).to_csv(Input_dir+'Exp.txt',header=False,index=False,sep='\t')
	pd.DataFrame(List).to_csv(Input_dir+'Symbol.txt', sep='\t', header=False, index=False)
	return List,A


def TF_expression(TFName,List,Input_dir,Match2,A):
    d= np.isin(TFName,List)
    TFName = TFName[d]
    List_idx=pd.DataFrame(range(len(List)),index=List)
    f=List_idx.loc[TFName][0].values
    TF = A[f, :]
    Match2 = Match2[np.isin(Match2[:, 1], TFName)]
    d = np.isin(TFName, Match2[:, 1])
    TFName = TFName[d]
    TF = TF[d, :]
    pd.DataFrame(TF).to_csv(Input_dir+'TFexp.txt', sep='\t', header=False, index=False)
    pd.DataFrame(TFName).to_csv(Input_dir+'TFName.txt', sep='\t', header=False, index=False)
    return TFName 

def index_generate(choosL_i,merged_s,merged_b,TFName):
    if choosL_i in merged_s.index:
        REid = merged_s.loc[choosL_i]['id_s']
        REid_b = merged_b.loc[choosL_i]['id_b']
    else:
        REid=''
        REid_b=''
    TFName_1 = np.delete(TFName, np.where(TFName == choosL_i))
    TFid = np.where(np.isin(TFName, TFName_1))[0]
    RE_s = '_'.join(map(str, REid))
    TF_s = '_'.join(map(str, TFid))
    RE_b = '_'.join(map(str, REid_b))
    return choosL_i, RE_s, TF_s, RE_b

def load_corr_RE_TG(Input_dir,List,Element_name,Element_name_bulk):
	Element_gene = pd.read_csv(Input_dir+"hg19_Peak_hg19_gene_u.txt", delimiter="\t", header=None)
	choosL = List
	index_ElementName=pd.DataFrame(np.arange(0, len(Element_name)),index=Element_name)
	index_Element_name_bulk=pd.DataFrame(np.arange(0, len(Element_name_bulk)),index=Element_name_bulk)
	index_Element_name_bulk=index_Element_name_bulk.groupby(index_Element_name_bulk.index).min()
	Element_gene.columns=['Element_name_b','Element_name_s','TG']
	Element_gene['value']=1
	#RE_all_s=Element_gene['Element_name_s'].unique()
	#RE_all_b=Element_gene['Element_name_b'].unique()
	#index_RE_all_s=pd.DataFrame(np.arange(0, len(RE_all_s)),index=RE_all_s)
	#index_RE_all_b=pd.DataFrame(np.arange(0, len(RE_all_b)),index=RE_all_b)
	Element_gene['id_s']=index_ElementName.loc[Element_gene['Element_name_s']][0].values
	Element_gene['id_b']=index_Element_name_bulk.loc[Element_gene['Element_name_b']][0].values
	merged_s = Element_gene.groupby('TG')['id_s'].agg(list).reset_index()
	merged_b = Element_gene.groupby('TG')['id_b'].agg(list).reset_index()
	merged_s = merged_s.set_index('TG')
	merged_b =merged_b.set_index('TG')
	#index_ElementName1=index_ElementName.loc[RE_all_s][0].values
	#index_Element_name_bulk1=index_Element_name_bulk.loc[RE_all_b][0].values
	return merged_s,merged_b

def load_motifbinding_chr(chrN,GRNdir,Input_dir,motifWeight):
    Motif_binding_temp=pd.read_csv(GRNdir+'MotifTarget_Matrix_'+chrN+'.txt',sep='\t',index_col=0)
    REs=Motif_binding_temp.index
    march_hg19_Regrion=pd.read_csv(Input_dir+'MotifTarget_hg19_hg38_'+chrN+'.txt',sep='\t',header=None)
    REoverlap=list(set(march_hg19_Regrion[1].values))
    Motif_binding_temp1=Motif_binding_temp.loc[REoverlap]
    REs=Motif_binding_temp1.index
    Motif_binding_temp=np.zeros([march_hg19_Regrion.shape[0],Motif_binding_temp.shape[1]])
    Motif_binding_temp=Motif_binding_temp1.loc[march_hg19_Regrion[1].values].values
    Motif_binding_temp=pd.DataFrame(Motif_binding_temp,index=march_hg19_Regrion[0].values,columns=Motif_binding_temp1.columns)
    Motif_binding_temp1=Motif_binding_temp.groupby(Motif_binding_temp.index).max()
    motifoverlap=list(set(Motif_binding_temp1.columns)&set(motifWeight.index))
    Motif_binding_temp1=Motif_binding_temp1[motifoverlap]
    motifWeight=motifWeight.loc[Motif_binding_temp1.columns]
    Motif_binding = np.diag(1.0 / (motifWeight.T + 0.1)) * Motif_binding_temp1.values.T
    Motif_binding = np.log1p(Motif_binding)
    return Motif_binding_temp1

def load_TFbinding(GRNdir,Input_dir,motifWeight,Match2,TFName,Element_name):
    from tqdm import tqdm
    motif_binding=pd.DataFrame()
    chrall=['chr'+str(i+1) for i in range(22)]
    chrall.append('chrX')
    for chrN in tqdm(chrall):
        Motif_binding_temp1=load_motifbinding_chr(chrN,GRNdir,Input_dir,motifWeight)
        motif_binding=pd.concat([motif_binding,Motif_binding_temp1],join='outer',axis=0)
    motif_binding=motif_binding.fillna(0)
    motif_binding=motif_binding.groupby(motif_binding.index).max()
    motifoverlap=list(set(motif_binding.columns)&set(motifWeight.index))
    Match2=Match2[np.isin(Match2[:, 0],motifoverlap), :]
    TF_binding_temp = np.zeros((len(TFName), len(Element_name)))
    Motif_binding=np.zeros((motif_binding.shape[1], len(Element_name)))
    Element_name_idx=pd.DataFrame(range(len(Element_name)),index=Element_name)
    idx=Element_name_idx.loc[motif_binding.index][0].values
    Motif_binding=np.zeros((motif_binding.shape[1], len(Element_name)))
    Motif_binding[:,idx]=motif_binding.loc[Element_name[idx]].values.T
    Motif_binding=pd.DataFrame(Motif_binding,index=motif_binding.columns,columns=Element_name)
    Match2=Match2[np.isin(Match2[:, 1],TFName), :]
    Motif_binding=Motif_binding.loc[Match2[:, 0]]
    Motif_binding.index=Match2[:, 1]
    TF_binding=Motif_binding.groupby(Motif_binding.index).sum()
    a = np.sum(TF_binding.values, axis=1)
    a[a == 0] =1
    TF_binding_n = np.diag(1.0 / a) @TF_binding.values
    TF_binding_n=pd.DataFrame(TF_binding_n.T,index=Element_name,columns=TF_binding.index)
    TF_binding=np.zeros((len(Element_name),len(TFName)))
    idx=np.isin(TFName,TF_binding_n.columns)
    TF_binding[:,idx]=TF_binding_n[TFName[idx]].values
    TF_binding=pd.DataFrame(TF_binding,index=Element_name,columns=TFName)
    TF_binding.to_csv(Input_dir+'TF_binding.txt',sep='\t',index=None,header=None)

def preprocess(RNA_file,ATAC_file,label_file,Input_dir,GRNdir,LINGER_dir,genome,method):
    if method=='LINGER':
        print('Overlap the regions with bulk data ...')
        subprocess.run(["sh", LINGER_dir+"extract_overlap_regions_LINGER.sh", Input_dir , GRNdir, genome])
        print('Generate psedubulk ...')
        TG_psedubulk,RE_psedubulk,indexfinal,clusternew=pseudo_bulk.pseudo_bulk(RNA_file,ATAC_file,label_file,Input_dir)
        print('Generate gene expression...')
        TFName = pd.read_csv(GRNdir+'TFName.txt',header=None)
        TFName.columns=['TFName']
        TFName=TFName['TFName'].values
        Match2=pd.read_csv(GRNdir+'Match2.txt',sep='\t')
        Match2=Match2.values
        List,A=gene_expression(Input_dir,GRNdir,TFName,Match2,TG_psedubulk)
        print('Generate TF expression...')
        TFName=TF_expression(TFName,List,Input_dir,Match2,A)
        print('Generate RE chromatin accessibility...')
        RE_psedubulk[RE_psedubulk > 100] = 100
        RE_psedubulk.to_csv(Input_dir+'Openness.txt',sep='\t',header=None,index=None)
        print('Generate TF binding...')
        Element_name_bulk = pd.read_csv(GRNdir+'all_hg19.txt', delimiter="\t", header=None)
        Element_name_bulk=Element_name_bulk[0].values
        Element_name = RE_psedubulk.index
        motifWeight=pd.read_csv(GRNdir+'motifWeight.txt',index_col=0,sep='\t')
        load_TFbinding(GRNdir,Input_dir,motifWeight,Match2,TFName,Element_name)
        print('Generate Index...')
    #Read hg19_Peak_hg19_gene_u.txt
        merged_s,merged_b=load_corr_RE_TG(Input_dir,List,Element_name,Element_name_bulk)
        from tqdm import tqdm
    #Assuming you have imported the necessary libraries and defined the variables
    #Create a progress bar for the loop
        choosL=List
        out=np.empty([len(choosL),4], dtype=object)
        for i in tqdm(range(len(choosL))):
            choosL_i = choosL[i]
            out[i, :] = index_generate(choosL_i,merged_s,merged_b,TFName)
        pd.DataFrame(out).to_csv(Input_dir+'index.txt', sep='\t', header=None, index=None)
    elif method=='baseline':
        print('Overlap the regions with bulk data ...')
        subprocess.run(["sh", LINGER_dir+"extract_overlap_regions_baseline.sh", Input_dir , GRNdir, genome])
    else:
        print('Method:' +method+ 'is not found! Please set method as baseline or LINGER')
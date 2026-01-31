import os
import numpy as np
import pandas as pd
#from LingerGRN.immupute_dis import immupute_dis
#import LingerGRN.pseudo_bulk as pseudo_bulk
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



def gene_expression(GRNdir,TG_pseudobulk,outdir):
	gene = pd.read_csv(GRNdir+'bulk_gene_all.txt')
	gene.columns=['gene']
	#gene=gene['gene']
	d1 = np.isin(TG_pseudobulk.index, gene['gene'].values)
	List = TG_pseudobulk.index[d1]
	A = np.log2(1 + TG_pseudobulk.loc[List])
	#Write Exp.txt and Symbol.txt
	pd.DataFrame(A).to_csv(outdir+'Exp.txt',sep='\t',index=False,header=False)
	pd.DataFrame(List).to_csv(outdir+'Symbol.txt', sep='\t', header=False, index=False)
	pd.DataFrame(A.columns).to_csv(outdir+'Col.txt',sep='\t',index=False,header=False)
	return List,A


def TF_expression(TFName,List,Match2,A,outdir):
    d= np.isin(TFName,List)
    TFName = TFName[d]
    List_idx=pd.DataFrame(range(len(List)),index=List)
    f=List_idx.loc[TFName][0].values
    TF = A.values[f, :]
    Match2 = Match2[np.isin(Match2[:, 1], TFName)]
    d = np.isin(TFName, Match2[:, 1])
    TFName = TFName[d]
    TF = TF[d, :]
    pd.DataFrame(TF).to_csv(outdir+'TFexp.txt', sep='\t', header=False, index=False)
    pd.DataFrame(TFName).to_csv(outdir+'TFName.txt', sep='\t', header=False, index=False)
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

def load_corr_RE_TG(List,Element_name,Element_name_bulk,outdir):
	Element_gene = pd.read_csv(outdir+"hg19_Peak_hg19_gene_u.txt", delimiter="\t", header=None)
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

def load_motifbinding_chr(chrN,GRNdir,motifWeight,outdir):
    Motif_binding_temp=pd.read_csv(GRNdir+'MotifTarget_Matrix_'+chrN+'.txt',sep='\t',index_col=0)
    REs=Motif_binding_temp.index
    march_hg19_Regrion=pd.read_csv(outdir+'MotifTarget_hg19_hg38_'+chrN+'.txt',sep='\t',header=None)
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

def load_TFbinding(GRNdir,motifWeight,Match2,TFName,Element_name,outdir):
    from tqdm import tqdm
    motif_binding=pd.DataFrame()
    chrall=['chr'+str(i+1) for i in range(22)]
    chrall.append('chrX')
    for chrN in tqdm(chrall):
        Motif_binding_temp1=load_motifbinding_chr(chrN,GRNdir,motifWeight,outdir)
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
    TF_binding.to_csv(outdir+'TF_binding.txt',sep='\t',index=None,header=None)

def extract_overlap_regions(genome,GRNdir,outdir,method):
    import pybedtools
    import pandas as pd
    import os
    os.makedirs(outdir, exist_ok=True)
    input_file = 'data/Peaks.txt'
    output_file = outdir+'Region.bed'
# Read the input file
    df = pd.read_csv(input_file, sep='\t',header=None)
    chromosomes = [item.split(':')[0] for item in df[0].values]
# Drop the first row
# Replace ':' and '-' with tabs
    df = df.replace({':': '\t', '-': '\t'}, regex=True)   
    chrall=['chr'+str(i+1) for i in range(23)]+['chrX']
    df=df[pd.DataFrame(chromosomes)[0].isin(chrall).values]
    df.to_csv(output_file, index=None, header=None)
    if method=='LINGER':
        if genome=='hg38':
            a = pybedtools.example_bedtool(outdir+'Region.bed')
            b = pybedtools.example_bedtool(GRNdir+'hg38_hg19_pair.bed')
            a_with_b = a.intersect(b, wa=True,wb=True)
            a_with_b.saveas(outdir+'temp.bed')
            a_with_b=pd.read_csv(outdir+'temp.bed',sep='\t',header=None)
            a_with_b[[6,7,8,0,1,2]].to_csv(outdir+'match_hg19_peak.bed',sep='\t',header=None,index=None)
        if genome=='hg19':
            a = pybedtools.example_bedtool(outdir+'Region.bed')
            b = pybedtools.example_bedtool(GRNdir+'hg19_hg38_pair.bed')
            a_with_b = a.intersect(b, wa=True,wb=True)
            a_with_b.saveas(outdir+'temp.bed')
            a_with_b=pd.read_csv(outdir+'temp.bed',sep='\t',header=None)
            a_with_b[[6,7,8,0,1,2]].to_csv(outdir+'match_hg19_peak.bed',sep='\t',header=None,index=None)
        a = pybedtools.example_bedtool(outdir+'match_hg19_peak.bed')
        b = pybedtools.example_bedtool(GRNdir+'RE_gene_corr_hg19.bed')
        a_with_b = a.intersect(b, wa=True,wb=True)
        a_with_b.saveas(outdir+'temp.bed')
        a_with_b=pd.read_csv(outdir+'temp.bed',sep='\t',header=None)
        a_with_b=a_with_b[(a_with_b[1].values==a_with_b[7].values)&(a_with_b[2].values==a_with_b[8].values)]
        a_with_b_n = pd.DataFrame({
        'column1': a_with_b[0] + ':' + a_with_b[1].astype(str) + '-' + a_with_b[2].astype(str),
        'column2': a_with_b[3] + ':' + a_with_b[4].astype(str) + '-' + a_with_b[5].astype(str),
        'column3': a_with_b[9]})
        a_with_b_n=a_with_b_n.drop_duplicates()
        a_with_b_n.to_csv(outdir+'hg19_Peak_hg19_gene_u.txt',sep='\t',header=None,index=None)
        chr_all=['chr'+str(i+1) for i in range(22)]
        chr_all.append('chrX')
        for chrtemp in chr_all:
            a = pybedtools.example_bedtool(outdir+'match_hg19_peak.bed')
            b = pybedtools.example_bedtool(GRNdir+'MotifTarget_matrix_'+chrtemp+'.bed')
            a_with_b = a.intersect(b, wa=True,wb=True)
            a_with_b.saveas(outdir+'temp.bed')
            a_with_b=pd.read_csv(outdir+'temp.bed',sep='\t',header=None)
            a_with_b=a_with_b[(a_with_b[1].values==a_with_b[7].values)&(a_with_b[2].values==a_with_b[8].values)]
            a_with_b_n = pd.DataFrame({
            'column1': a_with_b[3] + ':' + a_with_b[4].astype(str) + '-' + a_with_b[5].astype(str),
            'column2': a_with_b[6] + ':' + a_with_b[7].astype(str) + '-' + a_with_b[8].astype(str)})
            a_with_b_n=a_with_b_n.drop_duplicates()
            a_with_b_n.to_csv(outdir+'MotifTarget_hg19_hg38_'+chrtemp+'.txt',sep='\t',header=None,index=None)
            a = pybedtools.example_bedtool(GRNdir+genome+'_Peaks_'+chrtemp+'.bed')
            b = pybedtools.example_bedtool(outdir+'Region.bed')
            a_with_b = a.intersect(b, wa=True,wb=True)
            a_with_b.saveas(outdir+'Region_overlap_'+chrtemp+'.bed')
    if method=='baseline': 
        chr_all=['chr'+str(i+1) for i in range(22)]
        chr_all.append('chrX')
        for chrtemp in chr_all:
            a = pybedtools.example_bedtool(GRNdir+genome+'_Peaks_'+chrtemp+'.bed')
            b = pybedtools.example_bedtool(outdir+'Region.bed')
            a_with_b = a.intersect(b, wa=True,wb=True)
            a_with_b.saveas(outdir+'Region_overlap_'+chrtemp+'.bed')


def preprocess(TG_pseudobulk,RE_pseudobulk,GRNdir,genome,method,outdir):
    #package_dir = os.path.dirname(os.path.abspath(__file__))
    if method=='LINGER':
        extract_overlap_regions(genome,GRNdir,outdir,method)
        print('Mapping gene expression...')
        TFName = pd.read_csv(GRNdir+'TFName.txt',header=None)
        TFName.columns=['TFName']
        TFName=TFName['TFName'].values
        Match2=pd.read_csv(GRNdir+'Match2.txt',sep='\t')
        Match2=Match2.values
        List,A=gene_expression(GRNdir,TG_pseudobulk,outdir)
        print('Generate TF expression...')
        TFName=TF_expression(TFName,List,Match2,A,outdir)
        print('Generate RE chromatin accessibility...')
        RE_pseudobulk.to_csv(outdir+'Openness.txt',sep='\t',header=None,index=None)
        print('Generate TF binding...')
        Element_name_bulk = pd.read_csv(GRNdir+'all_hg19.txt', delimiter="\t", header=None)
        Element_name_bulk=Element_name_bulk[0].values
        Element_name = RE_pseudobulk.index
        motifWeight=pd.read_csv(GRNdir+'motifWeight.txt',index_col=0,sep='\t')
        load_TFbinding(GRNdir,motifWeight,Match2,TFName,Element_name,outdir)
        print('Generate Index...')
    #Read hg19_Peak_hg19_gene_u.txt
        merged_s,merged_b=load_corr_RE_TG(List,Element_name,Element_name_bulk,outdir)
        from tqdm import tqdm
    #Assuming you have imported the necessary libraries and defined the variables
    #Create a progress bar for the loop
        choosL=List
        out=np.empty([len(choosL),4], dtype=object)
        for i in tqdm(range(len(choosL))):
            choosL_i = choosL[i]
            out[i, :] = index_generate(choosL_i,merged_s,merged_b,TFName)
        pd.DataFrame(out).to_csv(outdir+'index.txt', sep='\t', header=None, index=None)
    elif method=='baseline':
        print('Overlap the regions with bulk data ...')
        #script_path = os.path.join( "extract_overlap_regions_baseline.sh")
        #subprocess.run(["sh", script_path, GRNdir, genome,outdir,workdir])
        extract_overlap_regions(genome,GRNdir,outdir,method)
    else:
        print('Method:' +method+ 'is not found! Please set method as baseline or LINGER')

import scanpy as sc
#set some figure parameters for nice display inside jupyternotebooks.
import scipy
import pandas as pd
import anndata
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import  csc_matrix
def get_adata(matrix,features,barcodes,label):
    ### generate the anndata
    matrix.data=matrix.data.astype(np.float32)
    adata=anndata.AnnData(X= csc_matrix(matrix.T))
    adata.var['gene_ids']=features[1].values
    adata.obs['barcode']=barcodes[0].values
    if len(barcodes[0].values[0].split("-"))==2:
        adata.obs['sample'] = [int(string.split("-")[1]) for string in barcodes[0].values]
    else:
        adata.obs['sample'] = 1
    rows_to_select=features[features[2]=='Gene Expression'].index
    adata_RNA = adata[:,rows_to_select]
    rows_to_select=features[features[2]=='Peaks'].index
    adata_ATAC = adata[:,rows_to_select]
### if you have the label (cell type annotation)
    idx=adata_RNA.obs['barcode'].isin(label['barcode_use'].values)
    adata_RNA=adata_RNA[idx]
    adata_ATAC=adata_ATAC[idx]
    label.index=label['barcode_use']
    adata_RNA.obs['label']=label.loc[adata_RNA.obs['barcode']]['label'].values
    #barcode_indices = np.where(np.isin(adata_RNA.obs['barcode'].values, label['barcode_use'].values))[0]
    #adata_ATAC = adata_ATAC[barcode_indices, :]
    adata_ATAC.obs['label']=label.loc[adata_ATAC.obs['barcode']]['label'].values
    adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
    adata_RNA, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < 5, :].copy()
    adata_RNA.var.index=adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids']=adata_RNA.var.index
    selected_barcode=list(set(adata_RNA.obs['barcode'].values)&set(adata_ATAC.obs['barcode'].values))
    barcode_idx=pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]
    barcode_idx=pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]
    return adata_RNA,adata_ATAC


def get_adata_h5(adata_RNA,adata_ATAC,label):
    ### generate the anndata
    if len(adata_RNA.obs['barcode'].values[0].split("-"))==2:
        adata_RNA.obs['sample'] = [int(string.split("-")[1]) for string in adata_RNA.obs['barcode'].values]
        adata_ATAC.obs['sample'] = [int(string.split("-")[1]) for string in adata_ATAC.obs['barcode'].values]
    else:
        adata_RNA.obs['sample'] = 1
        adata_ATAC.obs['sample'] = 1
### if you have the label (cell type annotation)
    idx=adata_RNA.obs['barcode'].isin(label['barcode_use'].values)
    adata_RNA=adata_RNA[idx]
    adata_ATAC=adata_ATAC[idx]
    label.index=label['barcode_use']
    adata_RNA.obs['label']=label.loc[adata_RNA.obs['barcode']]['label'].values
    #barcode_indices = np.where(np.isin(adata_RNA.obs['barcode'].values, label['barcode_use'].values))[0]
    #adata_ATAC = adata_ATAC[barcode_indices, :]
    adata_ATAC.obs['label']=label.loc[adata_ATAC.obs['barcode']]['label'].values
    adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
    adata_RNA, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < 5, :].copy()
    adata_RNA.var.index=adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids']=adata_RNA.var.index
    selected_barcode=list(set(adata_RNA.obs['barcode'].values)&set(adata_ATAC.obs['barcode'].values))
    barcode_idx=pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]
    barcode_idx=pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]
    return adata_RNA,adata_ATAC
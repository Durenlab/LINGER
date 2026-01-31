import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn import functional as F
from scipy.stats import pearsonr
from scipy.stats import spearmanr
#load data
import numpy as np
import pandas as pd
import random
from torch.optim import Adam
import os
from sklearn.linear_model import ElasticNet
from sklearn.datasets import make_regression
from sklearn.model_selection import KFold
import shap
hidden_size  = 64
hidden_size2 = 16
output_size = 1
from joblib import Parallel,delayed
seed_value = 42

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

#EWC
def EWC(fisher,params,net):
    params_n = list(net.parameters())
    EWC=0
    i=0 
    p=params_n[0]
    cost=(p-params[i])*fisher*(p-params[i])
    EWC=EWC+cost.sum()
    return EWC


def sc_nn(ii,gene_chr,TFindex,TFindex_bulk,REindex,REindex_bulk,REindex_bulk_match,Target,netall,adj_matrix_all,Exp,TF_match,input_size_all,fisherall,Opn,l1_lambda,fisher_w,activef):
    warnings.filterwarnings("ignore")
    alpha = 1
    eps=1e-12
    alpha = torch.tensor(alpha,dtype=torch.float32)
    gene_idx=gene_chr['id_s'].values[ii]-1
    gene_idx_b=int(gene_chr['id_b'].values[ii])-1
    TFidxtemp=TFindex[gene_idx]
    TFidxtemp=TFidxtemp.split('_')
    TFidxtemp=[int(TFidxtemp[k])+1 for k in range(len(TFidxtemp))]
    TFidxtemp_b=TFindex_bulk[gene_idx_b]
    TFidxtemp_b=TFidxtemp_b.split('_')
    TFidxtemp_b=[int(TFidxtemp_b[k]) for k in range(len(TFidxtemp_b))]
    TFtemp=Exp[np.array(TFidxtemp)-1,:]
    REidxtemp=REindex[gene_idx]
    REidxtemp_b_m=REindex_bulk_match[gene_idx]
    REidxtemp_b=REindex_bulk[gene_idx_b]
    REidxtemp=str(REidxtemp).split('_')
    REidxtemp_b_m=str(REidxtemp_b_m).split('_')
    REidxtemp_b=str(REidxtemp_b).split('_')
    if (len(REidxtemp)==1)&(REidxtemp[0]=='nan'):
        REidxtemp=[]
        REidxtemp_b_m=[]
        inputs=TFtemp+1-1
        L=np.zeros([len(TFidxtemp)+len(REidxtemp),len(TFidxtemp)+len(REidxtemp)])
        L=torch.tensor(L, dtype=torch.float32)
    else:
        REidxtemp=[int(REidxtemp[k])+1 for k in range(len(REidxtemp))]
        REidxtemp_b_m=[int(REidxtemp_b_m[k])+1 for k in range(len(REidxtemp_b_m))]
        REtemp=Opn[np.array(REidxtemp)-1,:]
        inputs=np.vstack((TFtemp, REtemp))
        adj_matrix=np.zeros([len(TFidxtemp)+len(REidxtemp),len(TFidxtemp)+len(REidxtemp)])
        AA=adj_matrix_all[np.array(REidxtemp)-1,:]
        AA=AA[:,np.array(TFidxtemp)-1]
        adj_matrix[:len(TFidxtemp),-len(REidxtemp):]=AA.T
        adj_matrix[-len(REidxtemp):,:len(TFidxtemp)]=AA
        A = torch.tensor(adj_matrix, dtype=torch.float32)
        D = torch.diag(A.sum(1))
        degree = A.sum(dim=1)
        degree += eps
        D_sqrt_inv = 1 / degree.sqrt()
        D_sqrt_inv = torch.diag(D_sqrt_inv)
        L = D_sqrt_inv@(D - A)@D_sqrt_inv
    if (len(REidxtemp_b)==1)&(REidxtemp_b[0]=='nan'):
        REidxtemp_b=[]
    else:
        REidxtemp_b=[int(REidxtemp_b[k]) for k in range(len(REidxtemp_b))]
    targets = torch.tensor(Target[gene_idx,:])
    inputs = torch.tensor(inputs,dtype=torch.float32)
    targets = targets.type(torch.float32)
    mean = inputs.mean(dim=1)
    std = inputs.std(dim=1)
    inputs = (inputs.T - mean) / (std+eps)
    inputs=inputs.T
    num_nodes=inputs.shape[0]
    y=targets.reshape(len(targets),1)     
    #trainData testData          
    input_size=int(input_size_all[gene_idx_b])
    loaded_net = Net(input_size,activef)
    loaded_net.load_state_dict(netall[gene_idx_b])
    params = list(loaded_net.parameters())
    fisher0=fisherall[gene_idx_b][0].data.clone()
    data0=pd.DataFrame(TFidxtemp)
    data1=pd.DataFrame(TFidxtemp_b)
    data0.columns=['TF']
    data1.columns=['TF']
    A=TF_match.loc[data0['TF'].values-1]['id_b']
    data0=pd.DataFrame(A)
    data0.columns=['TF']
    data1['id_b']=data1.index
    data0['id_s']=range(0,len(A))
    merge_TF=pd.merge(data0,data1,how='left',on='TF')
    if (len(REidxtemp)>0)&(len(REidxtemp_b)>0):
        data0=pd.DataFrame(REidxtemp_b_m)
        data1=pd.DataFrame(REidxtemp_b)
        data0.columns=['RE']
        data1.columns=['RE']
        data0['id_s']=data0.index
        data1['id_b']=data1.index
        merge_RE=pd.merge(data0,data1,how='left',on='RE')
        if merge_RE['id_b'].isna().sum()==0:
            good=1
            indexall=merge_TF['id_b'].values.tolist()+(merge_RE['id_b'].values+merge_TF.shape[0]).tolist()
        else: 
            good=0
    else:
        indexall=merge_TF['id_b'].values.tolist()
        good=1
    if good==1:      
        fisher=fisher0[:,np.array(indexall,dtype=int)]
        params_bulk = params[0][:,np.array(indexall,dtype=int)]
        with torch.no_grad():
            params_bulk = params_bulk.detach()     
        num_nodes=inputs.shape[0]
        n_folds = 5
        kf = KFold(n_splits=n_folds,shuffle=True,random_state=0)
        fold_size = len(inputs.T) // n_folds
        input_size = num_nodes
        mse_loss = nn.MSELoss()
        y_pred_all=0*(y+1-1)
        y_pred_all1=0*(y+1-1)
        y_pred_all1=y_pred_all1.numpy().reshape(-1)
        X_tr = inputs.T
        y_tr = y
        torch.manual_seed(seed_value)
        net = Net(input_size,activef)
        optimizer = Adam(net.parameters(),lr=0.01,weight_decay=l1_lambda)   
            #optimizer = Adam(net.parameters(),weight_decay=1)
            # Perform backpropagation
        Loss0=np.zeros([100,1])
        for i in range(100):
            # Perform forward pass
            y_pred = net(X_tr)
            # Calculate loss
            l1_norm = sum(torch.linalg.norm(p, 1) for p in net.parameters())
            #loss_EWC=EWC(fisher,params_bulk,net);
            l2_bulk = -1* fisher_w*  sum(sum(torch.mul(params_bulk,net.fc1.weight)))
            lap_reg = alpha * torch.trace(torch.mm(torch.mm(net.fc1.weight, L), net.fc1.weight.t()))
            loss = mse_loss(y_pred, y_tr) +l1_norm*l1_lambda+l2_bulk+lap_reg 
            Loss0[i,0]=loss.detach().numpy()
            # Perform backpropagation
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        np.random.seed(42)
        background = X_tr[np.random.choice(X_tr.shape[0], 50, replace=False)]
        explainer = shap.DeepExplainer(net,background)
        shap_values = explainer.shap_values(X_tr)
        warnings.resetwarnings()
        return net,shap_values,0.5,0.5,1,Loss0
    else:
        warnings.resetwarnings()
        return 0,0,0,0,0,0

def get_TSS(GRNdir,genome,TSS_dis):
    #import pyensembl
# Initialize Ensembl database for the desired genome assembly
    #ensembl = pyensembl.EnsemblRelease(release=release, species=species)  # For hg19
# ensembl = pyensembl.EnsemblRelease(release=104, species='mouse')  # For mm10
# Get all genes in the genome
    #genes = ensembl.genes()
# Retrieve TSS positions for each gene and store them in a list
    #tss_positions = []
    #strand=[]
    #chrom=[]
    #genesymbol=[]
    #for gene in genes:
        #tss_positions.append(gene.transcripts[0].start)
        #strand.append(gene.strand)
        #chrom.append('chr'+gene.contig)
        #genesymbol.append(gene.name)
    import pandas as pd
    Tssdf = pd.read_csv(GRNdir+'TSS_'+genome+'.txt',sep='\t',header=None)
    Tssdf.columns=['chr','TSS','symbol','strand']
    Tssdf['1M-']=Tssdf['TSS']-TSS_dis
    Tssdf['1M+']=Tssdf['TSS']+TSS_dis
    temp=Tssdf['1M-'].values
    temp[temp<1]=1
    Tssdf['1M-']=temp
    Tssdf=Tssdf[Tssdf['symbol']!='']
    Tssdf[['chr','1M-','1M+','symbol','TSS', 'strand']].to_csv('data/TSS_extend_1M.txt',sep='\t',index=None)

def load_data(GRNdir,outdir):
    gene_all=pd.DataFrame([])
    for i in range(22):
        chr='chr'+str(i+1)
        gene_file=GRNdir+chr+'_gene.txt'
        data0=pd.read_csv(gene_file,sep='\t',header=None)
        data0['chr']=chr
        data0['id_b']=data0.index+1
        gene_all=pd.concat([gene_all,data0])
    chr='chrX'
    gene_file=GRNdir+chr+'_gene.txt'
    data0=pd.read_csv(gene_file,sep='\t',header=None)
    data0['chr']=chr
    data0['id_b']=data0.index+1
    gene_all=pd.concat([gene_all,data0])
    gene_file=outdir+'Symbol.txt'
    data0=pd.read_csv(gene_file,sep='\t',header=None)
    data0.columns=['Symbol']
    data0['id_s']=data0.index+1
    gene_all.columns=['Symbol','chr','id_b']
    data_merge=pd.merge(data0,gene_all,how='left',on='Symbol')
    TFName_b=pd.read_csv(GRNdir+'TFName.txt',header=None,sep='\t')
    TFName_s=pd.read_csv(outdir+'TFName.txt',header=None,sep='\t')
    TFName_b.columns=['TF']
    TFName_s.columns=['TF']
    TFName_b['id_b']=TFName_b.index+1# index from 1
    TFName_s['id_s']=TFName_s.index+1# index from 1
    TF_match=pd.merge(TFName_s,TFName_b,how='left',on='TF')
    Opn_file=outdir+'Openness.txt'
    idx_file=outdir+'index.txt'
    geneexp_file=outdir+'Exp.txt'
    Target=pd.read_csv(geneexp_file,header=None,sep='\t')
    Target=Target.values
    #def sc_NN(gene_file,Opn_file,idx_file,geneexp_file,out_PCC,out_net):
    #alpha = torch.tensor(alpha,dtype=torch.float32)
    bind_file=outdir+'TF_binding.txt'
    adj_matrix_all=pd.read_csv(bind_file,header=None,sep='\t')
    adj_matrix_all=adj_matrix_all.values
    TFExp_file=outdir+'TFexp.txt'
    Opn=pd.read_csv(Opn_file,header=None,sep='\t')
    Opn=Opn.values
    idx=pd.read_csv(idx_file,header=None,sep='\t')
    Exp=pd.read_csv(TFExp_file,header=None,sep='\t')
    Exp=Exp.values
    return Exp,idx,Opn,adj_matrix_all,Target,data_merge,TF_match
def sc_nn_NN(ii,RE_TGlink_temp,Target,Exp,Opn,l1_lambda,activef):
    warnings.filterwarnings("ignore")
    alpha = 1
    eps=1e-12
    alpha = torch.tensor(alpha,dtype=torch.float32)
    if RE_TGlink_temp[0] in Exp.index:
        TFtemp = Exp.drop([RE_TGlink_temp[0]]).values
    else:
        TFtemp=Exp.values
    REtemp=Opn.loc[RE_TGlink_temp[1]].values
    inputs=np.vstack((TFtemp, REtemp))
    targets = torch.tensor(Target.loc[RE_TGlink_temp[0],:])
    inputs = torch.tensor(inputs,dtype=torch.float32)
    targets = targets.type(torch.float32)
    mean = inputs.mean(dim=1)
    std = inputs.std(dim=1)
    inputs = (inputs.T - mean) / (std+eps)
    inputs=inputs.T
    num_nodes=inputs.shape[0]
    y=targets.reshape(len(targets),1)     
    #trainData testData          
    input_size=int(num_nodes)
    mse_loss = nn.MSELoss()
    y_pred_all=0*(y+1-1)
    y_pred_all1=0*(y+1-1)
    y_pred_all1=y_pred_all1.numpy().reshape(-1)
    X_tr = inputs.T
    y_tr = y
    torch.manual_seed(seed_value)
    net = Net(input_size,activef)
    optimizer = Adam(net.parameters(),lr=0.01,weight_decay=l1_lambda)   
            #optimizer = Adam(net.parameters(),weight_decay=1)
            # Perform backpropagation
    Loss0=np.zeros([100,1])
    for i in range(100):
            # Perform forward pass
        y_pred = net(X_tr)
            # Calculate loss
        l1_norm = sum(torch.linalg.norm(p, 1) for p in net.parameters())
            #loss_EWC=EWC(fisher,params_bulk,net);
        #l2_bulk = -1* fisher_w*  sum(sum(torch.mul(params_bulk,net.fc1.weight)))
        #lap_reg = alpha * torch.trace(torch.mm(torch.mm(net.fc1.weight, L), net.fc1.weight.t()))
        loss = mse_loss(y_pred, y_tr) +l1_norm*l1_lambda#+l2_bulk+lap_reg 
        Loss0[i,0]=loss.detach().numpy()
            # Perform backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    np.random.seed(42)
    background = X_tr[np.random.choice(X_tr.shape[0], 50, replace=False)]
    explainer = shap.DeepExplainer(net,background)
    shap_values = explainer.shap_values(X_tr)
    warnings.resetwarnings()
    return net,shap_values,Loss0
    

def load_data_scNN(GRNdir,species):
    import pandas as pd
    if species=='New':
        Match2=pd.read_csv(GRNdir+'MotifMatch.txt',header=0,sep='\t')
    else:
        Match2=pd.read_csv(GRNdir+'Match_TF_motif_'+species+'.txt',header=None,sep='\t')
        Match2.columns = ['Motif','TF']
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

def RE_TG_dis(outdir):
    import pandas as pd
    import pybedtools
    import numpy as np
    print('Overlap the regions with gene loc ...')
    import os# Create the directory
    current_directory = os.getcwd()
    os.makedirs(outdir, exist_ok=True)
    import pandas as pd
    peakList=pd.read_csv(current_directory+'/data/Peaks.txt',index_col=None,header=None)
    peakList1=[temp.split(':')[0] for temp in peakList[0].values.tolist()]
    peakList2=[temp.split(':')[1].split('-')[0] for temp in peakList[0].values.tolist()]
    peakList3=[temp.split(':')[1].split('-')[1] for temp in peakList[0].values.tolist()]
    peakList['chr']=peakList1
    peakList['start']=peakList2
    peakList['end']=peakList3
    peakList[['chr','start','end']].to_csv(current_directory+'/data/Peaks.bed',sep='\t',header=None,index=None)
    TSS_1M=pd.read_csv(current_directory+'/data/TSS_extend_1M.txt',sep='\t',header=0)
    TSS_1M.to_csv(current_directory+'/data/TSS_extend_1M.bed',sep='\t',header=None,index=None)
    a = pybedtools.example_bedtool(current_directory+'/data/Peaks.bed')
    b = pybedtools.example_bedtool(current_directory+'/data/TSS_extend_1M.bed')
    a_with_b = a.intersect(b, wa=True,wb=True)
    a_with_b.saveas(outdir+'temp.bed')
    a_with_b=pd.read_csv(outdir+'temp.bed',sep='\t',header=None)
    a_with_b['RE']=a_with_b[0].astype(str) + ':' + a_with_b[1].astype(str) + '-' + a_with_b[2].astype(str)
    temp=a_with_b[['RE',6]]
    temp.columns=[['RE','gene']]
    temp['distance']=np.abs(a_with_b[7]-a_with_b[1])
    temp.to_csv(current_directory+'/data/RE_gene_distance.txt',sep='\t',index=None)
    
from tqdm import tqdm
import warnings
import time
import pandas as pd
import numpy as np
def training(GRNdir,method,outdir,activef,species):
    if method=='LINGER':
        hidden_size  = 64
        hidden_size2 = 16
        output_size = 1
        l1_lambda = 0.01 
        alpha_l = 0.01#elastic net parameter
        lambda0 = 0.00 #bulk
        fisher_w=0.1
        n_jobs=16

        Exp,idx,Opn,adj_matrix_all,Target,data_merge,TF_match=load_data(GRNdir,outdir)
        data_merge.to_csv(outdir+'data_merge.txt',sep='\t')
        chrall=['chr'+str(i+1) for i in range(22)]
        chrall.append('chrX')
        import warnings
        import time
        from tqdm import tqdm
        for i in range(23):
            netall_s={}
            shapall_s={}
            result=np.zeros([data_merge.shape[0],2])
            Lossall=np.zeros([data_merge.shape[0],100])
            chr=chrall[i]
            print(chr)
            idx_file1=GRNdir+chr+'_index.txt'
            idx_file_all=GRNdir+chr+'_index_all.txt'
            idx_bulk=pd.read_csv(idx_file1,header=None,sep='\t') 
            idxRE_all=pd.read_csv(idx_file_all,header=None,sep='\t')
            gene_chr=data_merge[data_merge['chr']==chr]
            N=len(gene_chr)
            TFindex=idx.values[:,2]
            REindex=idx.values[:,1]
            REindex_bulk_match=idx.values[:,3]
            REindex_bulk=idxRE_all.values[:,0]
            TFindex_bulk=idx_bulk.values[:,2]
            input_size_all=idx_bulk.values[:,3]
            fisherall = torch.load(GRNdir+'fisher_'+chr+'.pt')
            netall=torch.load(GRNdir+'all_models_'+chr+'.pt')
            
            for ii in tqdm(range(N)):
                warnings.filterwarnings("ignore")
                res=sc_nn(ii,gene_chr,TFindex,TFindex_bulk,REindex,REindex_bulk,REindex_bulk_match,Target,netall,adj_matrix_all,Exp,TF_match,input_size_all,fisherall,Opn,l1_lambda,fisher_w,activef)
                warnings.resetwarnings()
                index_all=gene_chr.index[ii]
                if res[4]==1:
                    result[index_all,0]=res[2]
                    result[index_all,1]=res[3]
                    netall_s[index_all]=res[0]
                    shapall_s[index_all]=res[1]
                    Lossall[index_all,:]=res[5].T
                else:
                    result[index_all,0]=-100
            result=pd.DataFrame(result)
            result.index=data_merge['Symbol'].values
            genetemp=data_merge[data_merge['chr']==chr]['Symbol'].values
            result=result.loc[genetemp]
            result.to_csv(outdir+'result_'+chr+'.txt',sep='\t')
            torch.save(netall_s,outdir+'net_'+chr+'.pt')
            torch.save(shapall_s,outdir+'shap_'+chr+'.pt')
            Lossall=pd.DataFrame(Lossall)
            Lossall.index=data_merge['Symbol'].values
            Lossall=Lossall.loc[genetemp]
            Lossall.to_csv(outdir+'Loss_'+chr+'.txt',sep='\t')            
    if method=='scNN':
        hidden_size  = 64
        hidden_size2 = 16
        output_size = 1
        l1_lambda = 0.01 
        alpha_l = 0.01#elastic net parameter
        lambda0 = 0.00 #bulk
        fisher_w=0.1
        n_jobs=16
        Exp,Opn,Target,RE_TGlink=load_data_scNN(GRNdir,species)
        import warnings
        import time
        from tqdm import tqdm
        netall_s={}
        shapall_s={}
        #result=np.zeros([data_merge.shape[0],2])
        chrall=[RE_TGlink[0][i][0].split(':')[0] for i in range(RE_TGlink.shape[0])]
        RE_TGlink['chr']=chrall
        chrlist=RE_TGlink['chr'].unique()
        for jj in tqdm(range(len(chrlist))):
            chrtemp=chrlist[jj]
            RE_TGlink1=RE_TGlink[RE_TGlink['chr']==chrtemp]
            Lossall=np.zeros([RE_TGlink1.shape[0],100])
            for ii in  range(RE_TGlink1.shape[0]):
                warnings.filterwarnings("ignore")
                #res = Parallel(n_jobs=n_jobs)(delayed(sc_nn_NN)(ii,RE_TGlink_temp,Target,netall,Exp,Opn,l1_lambda,activef)  for ii in tqdm(range(RE_TGlink.shape[0]))
                RE_TGlink_temp=RE_TGlink1.values[ii,:]
                res=sc_nn_NN(ii,RE_TGlink_temp,Target,Exp,Opn,l1_lambda,activef)
                warnings.resetwarnings()
                netall_s[ii]=res[0]
                shapall_s[ii]=res[1]
                Lossall[ii,:]=res[2].T   
            torch.save(netall_s,outdir+chrtemp+'_net.pt')
            torch.save(shapall_s,outdir+chrtemp+'_shap.pt')
            Lossall=pd.DataFrame(Lossall)
            Lossall.index=RE_TGlink1['gene'].values
            Lossall.to_csv(outdir+chrtemp+'_Loss.txt',sep='\t') 
        RE_TGlink.to_csv(outdir+'RE_TGlink.txt',sep='\t',index=None)


def get_TSS_ensembl(genome_short,gtf_file,GRNdir):
    import pyensembl
    import subprocess
    from pyensembl import Genome
    ensembl = Genome(
    reference_name=genome_short,
    annotation_name="My_annotation",
    gtf_path_or_url=gtf_file)
    ensembl.index()
    genes = ensembl.genes()
# Retrieve TSS positions for each gene and store them in a list
    tss_positions = []
    strand=[]
    chrom=[]
    genesymbol=[]
    for gene in genes:
        tss_positions.append(gene.transcripts[0].start)
        strand.append(gene.strand)
        chrom.append('chr'+gene.contig)
        genesymbol.append(gene.name)
    import pandas as pd
    Tssdf = pd.DataFrame({'chr': chrom, 'TSS': tss_positions, 'symbol': genesymbol,'strand': strand})
    Tssdf.to_csv(GRNdir+'TSS_'+genome_short+'.txt',sep='\t',index=None,header=0)
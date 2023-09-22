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

#EWC
def EWC(fisher,params,net):
    params_n = list(net.parameters())
    EWC=0
    i=0 
    p=params_n[0]
    cost=(p-params[i])*fisher*(p-params[i])
    EWC=EWC+cost.sum()
    return EWC


def sc_nn(ii,gene_chr,TFindex,TFindex_bulk,REindex,REindex_bulk,REindex_bulk_match,Target,netall,adj_matrix_all,Exp,TF_match,input_size_all,fisherall,Opn,l1_lambda,fisher_w):
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
    loaded_net = Net(input_size)
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
        net = Net(input_size)
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
        return net,shap_values,0.5,0.5,1,Loss0
    else:
        return 0,0,0,0,0,0


def load_data(GRNdir,Input_dir):
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
    gene_file=Input_dir+'Symbol.txt'
    data0=pd.read_csv(gene_file,sep='\t',header=None)
    data0.columns=['Symbol']
    data0['id_s']=data0.index+1
    gene_all.columns=['Symbol','chr','id_b']
    data_merge=pd.merge(data0,gene_all,how='left',on='Symbol')
    TFName_b=pd.read_csv(GRNdir+'TFName.txt',header=None,sep='\t')
    TFName_s=pd.read_csv(Input_dir+'TFName.txt',header=None,sep='\t')
    TFName_b.columns=['TF']
    TFName_s.columns=['TF']
    TFName_b['id_b']=TFName_b.index+1# index from 1
    TFName_s['id_s']=TFName_s.index+1# index from 1
    TF_match=pd.merge(TFName_s,TFName_b,how='left',on='TF')
    Opn_file=Input_dir+'Openness.txt'
    idx_file=Input_dir+'index.txt'
    geneexp_file=Input_dir+'Exp.txt'
    Target=pd.read_csv(geneexp_file,header=None,sep='\t')
    Target=Target.values
    #def sc_NN(gene_file,Opn_file,idx_file,geneexp_file,out_PCC,out_net):
    #alpha = torch.tensor(alpha,dtype=torch.float32)
    bind_file=Input_dir+'TF_binding.txt'
    adj_matrix_all=pd.read_csv(bind_file,header=None,sep='\t')
    adj_matrix_all=adj_matrix_all.values
    TFExp_file=Input_dir+'TFexp.txt'
    Opn=pd.read_csv(Opn_file,header=None,sep='\t')
    Opn=Opn.values
    idx=pd.read_csv(idx_file,header=None,sep='\t')
    Exp=pd.read_csv(TFExp_file,header=None,sep='\t')
    Exp=Exp.values
    return Exp,idx,Opn,adj_matrix_all,Target,data_merge,TF_match

from tqdm import tqdm
import warnings
import time
import pandas as pd
import numpy as np
def training(GRNdir,Input_dir,method):
    if method=='LINGER':
        hidden_size  = 64
        hidden_size2 = 16
        output_size = 1
        l1_lambda = 0.01 
        alpha_l = 0.01#elastic net parameter
        lambda0 = 0.00 #bulk
        fisher_w=0.1
        n_jobs=8
        netall_s={}
        shapall_s={}
        Exp,idx,Opn,adj_matrix_all,Target,data_merge,TF_match=load_data(GRNdir,Input_dir)
        data_merge.to_csv(Input_dir+'data_merge.txt',sep='\t')
        chrall=['chr'+str(i+1) for i in range(22)]
        chrall.append('chrX')
        import warnings
        import time
        from tqdm import tqdm
        for i in range(23):
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
                index_all=gene_chr.index[ii]
                warnings.filterwarnings("ignore")
                res=sc_nn(ii,gene_chr,TFindex,TFindex_bulk,REindex,REindex_bulk,REindex_bulk_match,Target,netall,adj_matrix_all,Exp,TF_match,input_size_all,fisherall,Opn,l1_lambda,fisher_w)
                warnings.resetwarnings()
                if res[4]==1:
                    result[index_all,0]=res[2]
                    result[index_all,1]=res[3]
                    netall_s[index_all]=res[0]
                    shapall_s[index_all]=res[1]
                    Lossall[index_all,]=res[5].T
                else:
                    result[index_all,0]=-100
            result=pd.DataFrame(result)
            result.index=data_merge['Symbol'].values
            genetemp=data_merge[data_merge['chr']==chr]['Symbol'].values
            result=result.loc[genetemp]
            result.to_csv(Input_dir+'result_'+chr+'.txt',sep='\t')
            torch.save(netall_s,Input_dir+'net_'+chr+'.pt')
            torch.save(shapall_s,Input_dir+'shap_'+chr+'.pt')
            Lossall=pd.DataFrame(Lossall)
            Lossall.index=data_merge['Symbol'].values
            Lossall=Lossall.loc[genetemp]
            Lossall.to_csv(Input_dir+'Loss_'+chr+'.txt',sep='\t') 


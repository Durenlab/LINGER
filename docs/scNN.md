# Other species tutorial
We support for the following species:
|genome_short  | species | species_ensembl | genome |
| --- | --- | --- | --- |
| canFam3 | dog | Canis_lupus_familiaris |CanFam3  | 
| danRer11|zebrafish|Danio_rerio|GRCz11|
|danRer10|zebrafish|Danio_rerio|GRCz10|
| dm6|fly|Drosophila_melanogaster|BDGP6|
| dm3|fly|Drosophila_melanogaster|BDGP5|
| rheMac8|rhesus|Macaca_mulatta|Mmul_8|
|mm10|mouse|Mus_musculus|GRCm38|
|mm9|mouse|Mus_musculus|NCBIM37|
|rn5|rat|Rattus_norvegicus|Rnor_5|
|rn6|rat|Rattus_norvegicus|Rnor_6|
|susScr3|pig|Sus_scrofa|Sscrofa10|
|susScr11|pig|Sus_scrofa|Sscrofa11|
|fr3|fugu|Takifugu_rubripes|FUGU5|
|xenTro9|frog|Xenopus_tropicalis|Xenopus_tropicalis_v9|
## Download the provided data 
We provide the TSS location for the above genome and the motif information.
```sh
Datadir=/path/to/LINGER/# the directory to store the data please use the absolute directory. Example: Datadir=/zfs/durenlab/palmetto/Kaya/SC_NET/code/github/combine/data/
mkdir $Datadir
cd $Datadir
wget --load-cookies /tmp/cookies.txt "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1V6Ds2P6SStLQJDpne-Ga-RkRZ9dfUjpR'  -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1V6Ds2P6SStLQJDpne-Ga-RkRZ9dfUjpR" -O provide_data.tar.gz && rm -rf /tmp/cookies.txt
```
or use the following link: [https://drive.google.com/file/d/1V6Ds2P6SStLQJDpne-Ga-RkRZ9dfUjpR/view?usp=sharing](https://drive.google.com/file/d/1V6Ds2P6SStLQJDpne-Ga-RkRZ9dfUjpR/view?usp=sharing)

Then unzip，
```sh
tar -xzf provide_data.tar.gz
```
## Prepare the input data
We take sc data of mm10 as an examle. The data is from the published paper (FOXA2 drives lineage plasticity and KIT pathway
activation in neuroendocrine prostate cancer).
The input data is the feature matrix from 10x sc-multiome data and Cell annotation/cell type label which includes: 
- Single-cell multiome data including matrix.mtx, features.tsv/features.txt, and barcodes.tsv/barcodes.txt
- Cell annotation/cell type label if you need the cell type-specific gene regulatory network (label.txt in our example).
<div style="text-align: right">
  <img src="barcode_mm10.png" alt="Image" width="300">
</div>  

If the input data is 10X h5 file or h5ad file from scanpy, please follow the instruction [h5/h5ad file as input](https://github.com/Durenlab/LINGER/blob/main/docs/h5_input.md) .

### sc data
We download the data using shell command line.
```sh
mkdir -p data
cd data
wget --load-cookies /tmp/cookies.txt "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1PDOmtO2oL-YVxKQY26jL91SAFedDALA0'  -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1PDOmtO2oL-YVxKQY26jL91SAFedDALA0" -O mm10_data.tar.gz && rm -rf /tmp/cookies.txt
tar -xzvf mm10_data.tar.gz
mv mm10_data/* ./
cd ../
```
We provide the cell annotation as following:
```sh
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=/1nFm5shjcDuDYhA8YGzAnYoYVQ_29_Yj4' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=/1nFm5shjcDuDYhA8YGzAnYoYVQ_29_Yj4" -O mm10_label.txt && rm -rf /tmp/cookies.txt
mv mm10_label.txt data/
```
## LINGER 
### Install
```sh
conda create -n LINGER python==3.10.0
conda activate LINGER
pip install LingerGRN==1.91
conda install bioconda::bedtools #Requirement
```
For the following step, we run the code in python.
#### Transfer the sc-multiome data to anndata  
We will transfer sc-multiome data to the anndata format and filter the cell barcode by the cell type label.
```python
import scanpy as sc
#set some figure parameters for nice display inside jupyternotebooks.
%matplotlib inline
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
#results_file = "scRNA/pbmc10k.h5ad"
import scipy
import pandas as pd
matrix=scipy.io.mmread('data/matrix.mtx')
features=pd.read_csv('data/features.txt',sep='\t',header=None)
barcodes=pd.read_csv('data/barcodes.txt',sep='\t',header=None)
label=pd.read_csv('data/mm10_label.txt',sep='\t',header=0)
from LingerGRN.preprocess import *
adata_RNA,adata_ATAC=get_adata(matrix,features,barcodes,label)# adata_RNA and adata_ATAC are scRNA and scATAC
```
#### Remove low counts cells and genes
```python
import scanpy as sc
sc.pp.filter_cells(adata_RNA, min_genes=200)
sc.pp.filter_genes(adata_RNA, min_cells=3)
sc.pp.filter_cells(adata_ATAC, min_genes=200)
sc.pp.filter_genes(adata_ATAC, min_cells=3)
selected_barcode=list(set(adata_RNA.obs['barcode'].values)&set(adata_ATAC.obs['barcode'].values))
barcode_idx=pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]
barcode_idx=pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]
```
#### Generate the pseudo-bulk/metacell:
```python
from LingerGRN.pseudo_bulk import *
samplelist=list(set(adata_ATAC.obs['sample'].values)) # sample is generated from cell barcode 
tempsample=samplelist[0]
TG_pseudobulk=pd.DataFrame([])
RE_pseudobulk=pd.DataFrame([])
singlepseudobulk = (adata_RNA.obs['sample'].unique().shape[0]*adata_RNA.obs['sample'].unique().shape[0]>100)
for tempsample in samplelist:
    adata_RNAtemp=adata_RNA[adata_RNA.obs['sample']==tempsample]
    adata_ATACtemp=adata_ATAC[adata_ATAC.obs['sample']==tempsample]
    TG_pseudobulk_temp,RE_pseudobulk_temp=pseudo_bulk(adata_RNAtemp,adata_ATACtemp,singlepseudobulk)                
    TG_pseudobulk=pd.concat([TG_pseudobulk, TG_pseudobulk_temp], axis=1)
    RE_pseudobulk=pd.concat([RE_pseudobulk, RE_pseudobulk_temp], axis=1)
    RE_pseudobulk[RE_pseudobulk > 100] = 100

import os
if not os.path.exists('data/'):
    os.mkdir('data/')
adata_ATAC.write('data/adata_ATAC.h5ad')
adata_RNA.write('data/adata_RNA.h5ad')
TG_pseudobulk=TG_pseudobulk.fillna(0)
RE_pseudobulk=RE_pseudobulk.fillna(0)
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv('data/Peaks.txt',header=None,index=None)
TG_pseudobulk.to_csv('data/TG_pseudobulk.tsv')
RE_pseudobulk.to_csv('data/RE_pseudobulk.tsv')
```
### Training model
Overlap the region with general GRN:
```python
Datadir='/path/to/LINGER/'# This directory should be the same as Datadir defined in the above 'Download the general gene regulatory network' section
GRNdir=Datadir+'provide_data/'
genome='mm10'
outdir='/path/to/output/' #output dir
activef='ReLU' 
method='scNN'
import torch
import subprocess
import os
import LingerGRN.LINGER_tr as LINGER_tr
LINGER_tr.get_TSS(GRNdir,genome,200000) # Here, 200000 represent the largest distance of regulatory element to the TG. Other distance is supported
LINGER_tr.RE_TG_dis(outdir)
```
Train for the LINGER model.
```python
import LingerGRN.LINGER_tr as LINGER_tr
activef='ReLU' # active function chose from 'ReLU','sigmoid','tanh'
LINGER_tr.training(GRNdir,method,outdir,activef)
```


### Cell population gene regulatory network
#### TF binding potential
The output is 'cell_population_TF_RE_binding.txt', a matrix of the TF-RE binding score.
```python
import LingerGRN.LL_net as LL_net
LL_net.TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
```

#### *cis*-regulatory network
The output is 'cell_population_cis_regulatory.txt' with 3 columns: region, target gene, cis-regulatory score.
```python
LL_net.cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
```
#### *trans*-regulatory network
The output is 'cell_population_trans_regulatory.txt', a matrix of the trans-regulatory score.
```python
LL_net.trans_reg(GRNdir,method,outdir,genome)
```

### Cell type sepecific gene regulaory network
There are 2 options:
1. infer GRN for a specific cell type, which is in the label.txt;
```python
celltype='1' #use a string to assign your cell type
```
2. infer GRNs for all cell types.
```python
celltype='all'
```
Please make sure that 'all' is not a cell type in your data.
#### Motif matching
Check whether homer is installed
```sh
which homer # run this in command line
```
If homer is not installed, use conda to install it
```sh
conda 
```
#### TF binding potential
The output is 'cell_population_TF_RE_binding_*celltype*.txt', a matrix of the TF-RE binding potential.
```python
LL_net.cell_type_specific_TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)# different from the previous version
```

#### *cis*-regulatory network
The output is 'cell_type_specific_cis_regulatory_{*celltype*}.txt' with 3 columns: region, target gene, cis-regulatory score.
```python
LL_net.cell_type_specific_cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir)
```

#### *trans*-regulatory network
The output is 'cell_type_specific_trans_regulatory_{*celltype*}.txt', a matrix of the trans-regulatory score.
```python
LL_net.cell_type_specific_trans_reg(GRNdir,adata_RNA,celltype,outdir)
```


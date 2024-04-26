# h5ad file as input
## case1. 10x filtered feature barcode matrix
### download h5 file
```sh
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```
### get the input data for LINGER
```python
import scanpy as sc
adata = scanpy.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', gex_only=False)
import scipy.sparse as sp
import pandas as pd
matrix=adata.X.T
adata.var['gene_ids']=adata.var.index
features=pd.DataFrame(adata.var['gene_ids'].values.tolist(),columns=[1])
features[2]=adata.var['feature_types'].values
barcodes=pd.DataFrame(adata.obs_names,columns=[0])
from LingerGRN.preprocess import *
adata_RNA,adata_ATAC=get_adata(matrix,features,barcodes,label)# adata_RNA and adata_ATAC are scRNA and scATAC
```
## case2. seperate RNA and ATAC h5ad file
### Read H5AD file as an AnnData object
```python
import scanpy as sc
adata_RNA = sc.read_h5ad('rna.h5ad')
adata_ATAC=sc.read_h5ad('ATAC.h5ad')
import pandas as pd
label=pd.read_csv('label.txt',sep='\t',header=0)
```
```python
adata_RNA
```

<div style="text-align: right">
  <img src="adata_RNA.png" alt="Image" width="500">
</div>

```python
adata_ATAC
```

<div style="text-align: right">
  <img src="adata_ATAC.png" alt="Image" width="500">
</div>

```python
label
```
<div style="text-align: right">
  <img src="label_PBMC.png" alt="Image" width="300">
</div>

### get the input data for LINGER

```python
import scipy.sparse as sp
matrix=sp.vstack([adata_RNA.X.T, adata_ATAC.X.T])
features=pd.DataFrame(adata_RNA.var['gene_ids'].values.tolist()+adata_ATAC.var['gene_ids'].values.tolist(),columns=[1])
K=adata_RNA.shape[1]
N=K+adata_ATAC.shape[1]
types = ['Gene Expression' if i <= K else 'Peaks' for i in range(0, N)]
features[2]=types
barcodes=pd.DataFrame(adata_RNA.obs['barcode'].values,columns=[0])
from LingerGRN.preprocess import *
adata_RNA,adata_ATAC=get_adata(matrix,features,barcodes,label)# adata_RNA and adata_ATAC are scRNA and scATAC
```
Then you could go back to the PBMC tutorial and continue with the 'Remove low counts cells and genes' step.

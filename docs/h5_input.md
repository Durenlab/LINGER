# h5ad file as input
## Read H5AD file as an AnnData object
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

```python
import scipy.sparse as sp
matrix=sp.vstack([adata_RNA.X.T, adata_ATAC.X.T])
features=pd.DataFrame(adata_RNA.var['gene_ids'].values.tolist()+adata_ATAC.var['gene_ids'].values.tolist(),columns=[1])
K=adata_RNA.shape[1]
N=K+adata_ATAC.shape[1]
types = ['Gene Expression' if i <= K else 'Peaks' for i in range(0, N)]
features[2]=types
barcodes=pd.DataFrame(adata_RNA.obs['barcode'].values,columns=[0])
```

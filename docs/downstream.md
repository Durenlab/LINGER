# Downstream analysis
## Regulatory Module
For this analysis, we first detect key TF-TG subnetworks (modules) from the cell population TF–TG trans-regulation. Then, we identify the differential regulatory modules differentially expressed from the case and control groups.
### Detect Module
#### Input
- pseudobulk gene expression: [TG_pseudobulk], 
- metadata including case and control in column 'group' and cell type annotation in column 'celltype': [metadata],
- LINGER outdir including a trans-regulatory network, 'cell_population_trans_regulatory.txt',
- GWAS data file, which is not necessary.

This is an example of the input.
```python
import pandas as pd
TG_pseudobulk = pd.read_csv('data/TG_pseudobulk.tsv',sep=',',header=0,index_col=0)
TG_pseudobulk = TG_pseudobulk[~TG_pseudobulk.index.str.startswith('MT-')] # remove the mitochondrion, if the specie is mouse, replace 'MT-' with 'mt-'
import scanpy as sc
adata_RNA = sc.read_h5ad('data/adata_RNA.h5ad')
label_all = adata_RNA.obs[['barcode','sample','label']]
label_all.index = label_all['barcode']
metadata = label_all.loc[TG_pseudobulk.columns]
metadata.columns = ['barcode','group','celltype']
outdir = 'output/'
GWASfile = ['AUD_gene.txt','AUD_gene2.txt']# GWAS file is a gene list with no head (Optional)
```
```python
TG_pseudobulk
```
```python
metadata
```
## Driver Score

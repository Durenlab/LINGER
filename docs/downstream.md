# Downstream analysis
## Regulatory Module
For this analysis, we first detect key TF-TG subnetworks (modules) from the cell population TF–TG trans-regulation. Then, we identify the differential regulatory modules differentially expressed from the case and control groups.
### Detect Module
#### Input
- pseudobulk gene expression: [TG_pseudobulk], please make sure the data is after removing batch effect
- metadata including case and control in column 'group' and cell type annotation in column 'celltype': [metadata]. Note that the case is 1 and the control is 1.
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
<div style="text-align: right">
  <img src="RNA_ds.jpg" alt="Image" width="600">
</div>

```python
metadata 
```
<div style="text-align: right">
  <img src="metadata_ds.jpg" alt="Image" width="300">
</div>

```python
K=10 #k is the number of modules, a tunning parameter
Module_result=Module_trans(outdir,metadata,TG_pseudobulk,K,GWASfile)
```
The output is Module_result object. There are 3 items in this object: 
- S_TG, which representing the module assigned for eahc gene;
- pvalue_all, the p-value of the differential module t-test comparing the case and control groups;
- t_value, the t-value of the t-test, positive value representing group 1 is more active, and negative value representing group 0 is more active.
```python
Module_result.S_TG
```
<div style="text-align: right">
  <img src="S_TG.png" alt="Image" width="100">
</div>

```python
Module_result.pvalue_all
```
<div style="text-align: right">
  <img src="pvalue_all.png" alt="Image" width="400">
</div>

```python
Module_result.tvalue_all
```
<div style="text-align: right">
  <img src="tvalue_all.png" alt="Image" width="400">
</div>
## Driver Score

# Downstream analysis
## Regulatory Module
For this analysis, we first detect key TF-TG subnetworks (modules) from the cell population TF–TG trans-regulation. Then, we identify the differential regulatory modules differentially expressed from the case and control samples/cell types.
### Detect Module
#### Input
- pseudobulk gene expression [TG_pseudobulk], 
- metadata including case and control['group'], and cell type annotation['celltype'] [metadata],
- LINGER outdir including trans regulatory network,
- GWAS data, which is not neccessary.
## Driver Score

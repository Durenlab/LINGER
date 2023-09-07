# Construct the gene regulatory netwok by intergrating the general network (no training)
## Instruction
This tutorial delineates an approach for constructing cell type-specific gene regulatory networks from single-cell data using feature engineering methodology. Rather than directly training models on the single cell data, this workflow offers a rapid approach.

Just as the following figure, we combine the single cell data ($O, E$, and $C$ in the figure) and the prior gene regulatory network structure with the parameter $\alpha,\beta,d,B$, and $\gamma$.
![Example Figure](https://github.com/Durenlab/LINGER/blob/main/feature_engineering.jpg)
In this tutorial, we will 1. load the data we provide, 2. preprocess, 3. prepare the input data. 4. generate the cell population level gene regulatory network, 5. generate the cell type specific gene regulatory network.
## Download the the general gene regulatory network 
We provide the general gene regulatory network
```
LINGERdir=/path/to/LINGER/
cd $LINGERdir
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1vM8btN3LWu699YiPH0JyjIV_LZBIxiX_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1vM8btN3LWu699YiPH0JyjIV_LZBIxiX_" -O data_bulk.tar.gz && rm -rf /tmp/cookies.txt
```
Then unzip，
```sh
tar -xzf data_bulk.tar.gz
```
## Prepare the input data
The input data should be in same directory (Input_dis in this tutorial), which includes: 
- Single-cell multiome data including gene expression (RNA.txt in our example) and chromatin accessibility (ATAC.txt in our example).
- Cell annotation/cell type label if you need the cell type specific gene regulatory network (label.txt in our example).
### RNA-seq
The row of RNA-seq is gene symbol; the column is barcode; the value is the count matrix. Here is our example:
![Image Alt Text](RNA.png)
### ATAC-seq
The row is regulatory element/genomic region; the column is barcode, which is the same order with RNA-seq data; the value is the count matrix. Here is our example:
![Image Alt Text](ATAC.png)
### Cell annotation/cell type label
The row is cell barcode, which is the same order with RNA-seq data; there is one column 'Annotation', which is the cell type label. It could be a number or the string. Here is our example:
![Image Alt Text](label.png)
## Gene regulatory network inference
### Preprocess
Map the regions to the given regions by running the following code in linux. The output is overlaped region in each chromtin (Region_overlap_chr*.bed) in the same directory of input data.
```sh
Input_dir=/path/to/dir/ # all the input file should be in this directory
genome=hg38 # only hg38 and hg19 supported
cd $Input_dir
cat ATAC.txt|cut -f 1 |sed '1d' |sed 's/:/\t/g'| sed 's/-/\t/g' > Region.bed
GRNdir=$LINGERdir/data_bulk
$LINGERdir/extract_overlap_regions.sh "$GRNdir" $genome
cd $LINGERdir
```
### Load the input data.
We load the input data and import the function. The following sections are in python.
```python
import LL_net
Input_dir='/zfs/durenlab/palmetto/Kaya/SC_NET/code/github/version1/Input/'
GRNdir='/zfs/durenlab/palmetto/Kaya/SC_NET/code/github/version1/data_bulk/'
RNA_file='RNA.txt'
labels='label.txt'
ATAC_file='ATAC.txt'
genome='hg38'
```
### cell population gene regulatory network
#### TF binding potential
The output is 'cell_population_TF_RE_binding.txt', a matrix of the TF-RE binding score.
```python
result=LL_net.TF_RE_binding(Input_dir,GRNdir,RNA_file,ATAC_file,genome)
result.to_csv(Input_dir+'cell_population_TF_RE_binding.txt',sep='\t')
```
#### *cis*-regulatory network
The output is 'cell_population_cis_regulatory.txt' with 3 columns: region, target gene, cis-regulatory score.
```python
cis=LL_net.cis_reg(Input_dir,GRNdir,RNA_file,ATAC_file,genome)
cis.to_csv(Input_dir+'cell_population_cis_regulatory.txt',sep='\t',header=None,index=None)
```
#### *trans*-regulatory network
The output is 'cell_population_trans_regulatory.txt', a matrix of the trans-regulatory score.
```python
trans=LL_net.trans_reg(Input_dir,GRNdir,RNA_file,labels,ATAC_file)
trans.to_csv(Input_dir+'cell_population_trans_regulatory.txt',sep='\t')
```
### cell type sepecific gene regulaory network
There are 2 options:
1. infer GRN for a specific cell type, which is in the label.txt;
```python
celltype='0'#use a string to assign your cell type
```
2. infer GRNs for all cell types.
```python
celltype='all'
```
Please make sure that 'all' is not a cell type in your data.
#### TF binding potential
The output is 'cell_population_TF_RE_binding_*celltype*.txt', a matrix of the TF-RE binding potential.
```python
LL_net.cell_type_specific_TF_RE_binding_chr(RNA_file,ATAC_file,labels,Input_dir,GRNdir,chrN,genome,celltype)
```
#### *cis*-regulatory network
The output is 'cell_type_specific_cis_regulatory_{*celltype*}.txt' with 3 columns: region, target gene, cis-regulatory score.
```python
cis=LL_net.cell_type_specific_cis_reg(Input_dir,GRNdir,RNA_file,ATAC_file,genome,celltype)
```
#### *trans*-regulatory network
The output is 'cell_type_specific_trans_regulatory_{*celltype*}.txt', a matrix of the trans-regulatory score.
```python
trans=LL_net.cell_type_specific_trans_reg(Input_dir,GRNdir,RNA_file,labels,ATAC_file,celltype)
```


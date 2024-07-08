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
wget --load-cookies /tmp/cookies.txt "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1lAlzjU5BYbpbr4RHMlAGDOh9KWdCMQpS'  -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lAlzjU5BYbpbr4RHMlAGDOh9KWdCMQpS" -O data_bulk.tar.gz && rm -rf /tmp/cookies.txt
```
or use the following link: [https://drive.google.com/file/d/1lAlzjU5BYbpbr4RHMlAGDOh9KWdCMQpS/view?usp=sharing](https://drive.google.com/file/d/1lAlzjU5BYbpbr4RHMlAGDOh9KWdCMQpS/view?usp=sharing)

Then unzip，
```sh
tar -xzf data_bulk.tar.gz
```
## Prepare the input data
We take sc data of mm10 as an examle. The data is from the published paper [].
The input data is the feature matrix from 10x sc-multiome data and Cell annotation/cell type label which includes: 
- Single-cell multiome data including matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz.
- Cell annotation/cell type label if you need the cell type-specific gene regulatory network (PBMC_label.txt in our example).
<div style="text-align: right">
  <img src="label_PBMC.png" alt="Image" width="300">
</div>  

If the input data is 10X h5 file or h5ad file from scanpy, please follow the instruction [h5/h5ad file as input](https://github.com/Durenlab/LINGER/blob/main/docs/h5_input.md) .



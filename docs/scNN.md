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


## Download the general gene regulatory network 
We provide the general gene regulatory network, please download the data first.
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

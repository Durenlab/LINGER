# LINGER
## Introduction
LINGER (LIfelong neural Network for GEne Regulation) is a novel method to infer GRNs from single-cell multiome data built on top of [PyTorch](https://pytorch.org/).

LINGER incorporates both 1) atlas-scale external bulk data across diverse cellular contexts and 2) the knowledge of transcription factor (TF) motif matching to cis-regulatory elements as a manifold regularization to address the challenge of limited data and extensive parameter space in GRN inference.
## Analysis tasks for single cell multiome data
- Infer gene regulatory network
- Benchmark gene regulatory network
- Explainable dimensionality reduction (transcription factor activity, availiable for single cell or bulk RNA-seq data)
- In silico pertubation

In the user guide, we provide an overview of each task. 
## Basic installation
LINGER can be installed by pip
```sh
conda create -n LINGER python==3.10.0
conda activate LINGER
pip install LingerGRN==1.52
conda install bioconda::bedtools # Requirment
```
## Documentation

We provide several tutorials and user guide. If you find our tool useful for your research, please consider citing the LINGER manuscript.

|                           |                           |                           |
|:-------------------------:|:-------------------------:|:-------------------------:|
| [User guide](https://github.com/Durenlab/LINGER/blob/main/docs/User_guide.md) | [PBMCs tutorial](https://github.com/Durenlab/LINGER/blob/main/docs/PBMC.md) |[H1 cell line tutorial](https://github.com/Durenlab/LINGER/blob/main/docs/GRN_infer.md)  |
|[GRN benchmark](https://github.com/Durenlab/LINGER/blob/main/docs/Benchmark.md)  | [In silico perturbation](https://github.com/Durenlab/LINGER/blob/main/docs/perturb.md) | [Other species](https://github.com/Durenlab/LINGER/blob/main/docs/scNN.md) |
    

## Reference
> If you use LINGER, please cite:
> 
> [Yuan, Qiuyue, and Zhana Duren. "Continuous lifelong learning for modeling of gene regulation from single cell multiome data by leveraging atlas-scale external data." bioRxiv (2023).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418251/)

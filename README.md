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
pip install LingerGRN==1.45
conda install bioconda::bedtools # Requirment
```
## Documentation

We provide several tutorials and user guide. If you find our tool useful for your research, please consider citing the LINGER manuscript.

|                           |                           |                           |
|:-------------------------:|:-------------------------:|:-------------------------:|
| [User guide](https://github.com/Durenlab/LINGER/blob/main/docs/User_guide.md) | [Construct the gene regulatory network](https://github.com/Durenlab/LINGER/blob/main/docs/PBMC.md) | [GRN benchmark](https://github.com/Durenlab/LINGER/blob/main/docs/Benchmark.md) |
| [Identify driver regulators by TF activity](https://github.com/Durenlab/LINGER/blob/main/docs/TFactivity.md) | [In silico perturbation](https://github.com/Durenlab/LINGER/blob/main/docs/perturb.md) | [PBMCs data turtorial](https://example.com) |
    

## Reference
> If you use LINGER, please cite:
> 
> [Yuan, Qiuyue, and Zhana Duren. "Continuous lifelong learning for modeling of gene regulation from single cell multiome data by leveraging atlas-scale external data." bioRxiv (2023).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418251/)

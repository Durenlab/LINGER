# LINGER
## Introduction
LINGER (LIfelong neural Network for GEne Regulation)is a novel method to infer GRNs from single-cell multiome data built on top of [PyTorch](https://pytorch.org/).

LINGER incorporates both 1) atlas-scale external bulk data across diverse cellular contexts and 2) the knowledge of transcription factor (TF) motif matching to cis-regulatory elements as a manifold regularization to address the challenge of limited data and extensive parameter space in GRN inference.
## Analysis tasks for single cell multiome data
- Infer gene regulatory network
- Benchmark gene regulatory network
- Explainable dimensionality reduction (transcription factor activity, availiable for single cell or bulk RNA-seq data)
- Explain GWAS SNPs
- In silico pertubation

In the user guide, we provide an overview of each task. 
## Basic installation
LINGER can be installed by pip
```sh
pip install LINGER==1.0
```
## Documentation

We provide several tutorials and user guide. If you find it useful for your research, please consider citing the LINGER manuscript.

|                           |                           |                           |
|:-------------------------:|:-------------------------:|:-------------------------:|
| [User guide](https://example.com) | [Construct the gene regulatory network](https://github.com/Durenlab/LINGER/blob/main/tutorial1.md) | [GRN benchmark](https://example.com) |
| [Identify driver regulators by TF activity](https://github.com/Durenlab/LINGER/blob/main/tutorial2.md) | [GWAS trait annotation](https://example.com) | [In silico perturbation](https://example.com) |
    

## Reference
> If you use LINGER, please cite:
> 
> [Yuan, Qiuyue, and Zhana Duren. "Continuous lifelong learning for modeling of gene regulation from single cell multiome data by leveraging atlas-scale external data." bioRxiv (2023).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10418251/)

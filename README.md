# LINGER
## Introduction
LINGER (LIfelong neural Network for GEne Regulation)is a novel deep learning-based method to infer GRNs from single-cell multiome data with paired gene expression and chromatin accessibility data from the same cell. LINGER incorporates both 1) atlas-scale external bulk data across diverse cellular contexts and 2) the knowledge of transcription factor (TF) motif matching to cis-regulatory elements as a manifold regularization to address the challenge of limited data and extensive parameter space in GRN inference.
## Requirements
We use the following Python packages with python 3.10.9: 
torch: 1.13.1+cu117; numpy: 1.23.5; scipy: 1.10.1; pandas: 1.5.3; sklearn: 1.2.1; shap: 0.41.0; joblib: 1.2.0; 
### Run
#### step 0, install the package
For the first step, we can download the code by
```sh
git clone https://github.com/Durenlab/LINGER.git
cd LINGER
```
Then download the datasets:
```sh
wget https://drive.google.com/file/d/1miQkV1mUjBa7wFoPKcKXQwHR9ReCCBfO/view?usp=sharing
wget https://drive.google.com/file/d/1Rj7RbzY-8Tc8sRWJ_dfbsUlWTx9PFX1A/view?usp=sharing
```


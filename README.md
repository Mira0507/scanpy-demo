# scanpy-demo

This repository contains demo scripts for performing single-cell RNA sequencing (scRNA-seq) 
analysis using [Scanpy](https://scanpy.readthedocs.io/en/stable/), a scalable Python-based 
framework for analyzing single-cell gene expression data.

## Repository Structure

```
$ tree
.
├── env.yaml
├── README.md
└── scanpy
    ├── helpers.R
    └── scripts
        ├── sampletable.tsv
        ├── scanpy-integrated-markers-all.Rmd
        ├── scanpy-integration-all.Rmd
        ├── scanpy-qc-unintegrated-all.Rmd
        └── scanpy-unintegrated-markers-all.Rmd
```

- `env.yaml`: conda environment
- `scanpy/helpers.R`: helper functions
- `scanpy/scripts/sampletable.tsv`: sampletable pointing to input `h5` files generated 
by 10x Genomics [CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/latest).
With no further modification, the parental directory of `outs/filtered_feature_bc_matrix.h5`
is pointed.
- `scanpy/scripts/scanpy-qc-unintegrated-all.Rmd`: 
    - Quality control (removing doublet and outlier cells)
    - Normalization
    - Dimensionality reduction
    - Clustering
- `scanpy/scripts/scanpy-unintegrated-markers-all.Rmd`: 
    - Marker gene computation
- `scanpy/scripts/scanpy-integration-all.Rmd`: 
    - Integration based on deep learning using
    [scVI](https://pubmed.ncbi.nlm.nih.gov/30504886/)
- `scanpy/scripts/scanpy-integrated-markers-all.Rmd`: 
    - Marker gene computation

## Setup 

1. Clone the repository

Clone the repository to your local working directory. If your authentication 
method for GitHub is [SSH (Secure Shell Protocol)](https://www.ssh.com/academy/ssh-keys), 
run the following:

```
git clone git@github.com:Mira0507/scanpy-demo.git 
cd scanpy-demo
```

Otherwise, clone using the web URL below:

```
git clone https://github.com/Mira0507/scanpy-demo.git
cd scanpy-demo
```

2. Create Conda environment

The package management of the current workflow relies on conda and mamba. 
Ensure you have conda and mamba ready in your terminal. For more information, 
refer to the following pages:

- [Conda documentation](https://docs.conda.io/projects/conda/en/stable/)
- [Mamba user guide](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html)

If you are ready to set up your main conda environment, follow the command below:

```
# Assume you are in the scanpy-demo directory
$ mamba env create --prefix ./env --file env.yaml
```

This will create a new conda environment named `env` in the current directory.


## Analysis


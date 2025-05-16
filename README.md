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
    - Integration based on variational autoencoder-based using
    [scVI](https://pubmed.ncbi.nlm.nih.gov/30504886/)
- `scanpy/scripts/scanpy-integrated-markers-all.Rmd`: 
    - Marker gene computation

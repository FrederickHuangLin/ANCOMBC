# ANCOMBC
R package source code for implementing Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC).

ANCOM-BC is a methodology of differential abundance (DA) analysis that is designed to determine taxa that are differentially abundant with respect to the covariate of interest. For more details, please refer to: https://www.nature.com/articles/s41467-020-17041-7

Currently, the ANCOMBC package has been accepted on Bioconductor and added to its nightly builds. This version of ANCOMBC is able to perform covariates adjustment and global tests as mentioned in this paper.

## To install the latest release version of ANCOMBC

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
```

## Usage

```r
library(ANCOMBC)
?ancombc 
```
Author: Huang Lin, Shyamal Das Peddada

Maintainer: Huang Lin: <HUL40@pitt.edu>

# ANCOMBC
R package source code for implementing Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC).

ANCOM-BC is a methodology of differential abundance (DA) analysis that is designed to determine taxa that are differentially abundant with respect to the covariate of interest. For more details, please refer to: https://www.nature.com/articles/s41467-020-17041-7

Currently, the ANCOMBC package has been accepted on Bioconductor and added to its nightly builds. This version of ANCOMBC is able to perform covariates adjustment and global tests as mentioned in this paper.

## To install the latest release version of ANCOMBC

1. In order to use the ‘bioc-devel’ version of *Bioconductor* during the mid-April to mid-October release cycle, use the release version of *R* and invoke the function `install(version="devel")`
2. In order to use the ‘bioc-devel’ version of *Bioconductor* during the mid-October to mid-April release cycle, you must install the devel version of *R* and invoke the function `install(version="devel")`
3. For details, see: https://bioconductor.org/developers/how-to/useDevel/

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("ANCOMBC")
```

## Usage

```r
library(ANCOMBC)
?ancombc 
```
Author: Huang Lin, Shyamal Das Peddada

Maintainer: Huang Lin: <HUL40@pitt.edu>

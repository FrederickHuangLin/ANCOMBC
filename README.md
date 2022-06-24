# ANCOMBC
ANCOMBC is a package for differential abundance ([ANCOM-BC](https://www.nature.com/articles/s41467-020-17041-7) and [ANCOM](https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663)) and correlation (SECOM (in review)) analyses for microbiome data. Microbiome data are typically subject to two sources of biases: unequal sampling fractions (sample-specific biases) and differential sequencing efficiencies (taxon-specific biases). ANCOMBC package includes methodologies that aim to correct these biases and construct statistically consistent estimators.

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
?ancom
?secom_linear
?secom_dist
```
Author: Huang Lin

Maintainer: Huang Lin: <huanglinfrederick@gmail.com>

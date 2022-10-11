# ANCOMBC

ANCOMBC is a package containing differential abundance (DA) and correlation 
analyses for microbiome data. Specifically, the package includes 
Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2, manuscript in preparation),
Analysis of Compositions of Microbiomes with Bias Correction ([ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7)), and 
Analysis of Composition of Microbiomes ([ANCOM](https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663)) for DA analysis, and Sparse 
Estimation of Correlations among Microbiomes ([SECOM](https://doi.org/10.1038/s41467-022-32243-x)) for correlation 
analysis. Microbiome data are typically subject to two sources of biases: 
unequal sampling fractions (sample-specific biases) and differential 
sequencing efficiencies (taxon-specific biases). Methodologies included in 
the ANCOMBC package are designed to correct these biases and construct 
statistically consistent estimators.

## To install the latest release version of ANCOMBC

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
```

## Usage

```r
library(ANCOMBC)
?ancombc2
?ancombc 
?ancom
?secom_linear
?secom_dist
```
Author: Huang Lin

Maintainer: Huang Lin: <huanglinfrederick@gmail.com>

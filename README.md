# ANCOMBC

*Author & Maintainer: Huang Lin: <huanglinfrederick@gmail.com>*

ANCOMBC is a package containing differential abundance (DA) and correlation 
analyses for microbiome data. Specifically, the package includes 
Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2, submitted),
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
## Commonly asked questions

**1. Q: What exactly are the differences between `formula` and `group` arguments?**

A: When using the `formula` argument in `ancombc` and `ancombc2`, it is important to include all variables of your experiment that can have an influence on the microbial abundances. The `group` argument is optional, and should only be included if you are interested in detecting structural zeros (presence/absence test) or performing multi-group comparisons, such as the global test in `ancombc` and `ancombc2`, as well as the pairwise directional test, Dunnett’s type of test, and trend test in `ancombc2`. For example, if your variable of interest is a continuous variable, such as `age`, and you have other categorical variables that need to be adjusted but are not your research interest, you can leave `group = NULL`. However, if your variable of interest is a categorical variable with more than three levels and you are interested in performing multi-group comparisons, you should specify it in both the `formula` and `group` arguments. The formula tells `ancombc` and `ancombc2` to perform bias-correction for the specified variables and generate the primary results without multi-group comparisons, while the `group` argument is used to conduct multi-group comparisons and correct p-values for multiple comparisons.

**2. Q: Why are some taxa absent from primary results?**

A: Firstly, in the analysis, taxa with prevalences below `prv_cut` will be excluded. In addition, if the taxa contain structural zeros, they will be considered significant only by the presence/absence test, and not by the ANCOM-BC or ANCOM-BC2 methodology. Therefore, we have decided to summarize these results separately and not include them in the primary results of `ancombc` or `ancombc2`. You can access the results of the presence/absence test by checking `zero_ind`.

**3. Q: Say, I have a `group` variable containing `A`, `B`, and `C`, what do `lfc_(Intercept)`, `lfc_groupB`, and `lfc_groupC`mean in the primary results?** 

A: `lfc_groupB`, and `lfc_groupC` are the log (natural log) fold-changes with respect to the reference group, which is group `A` by default (you can use the `factor` function in R to change the reference group). Thus, they mean `LFC (group B - group A)` and `LFC (group C - group A)`,respectively. `lfc_(Intercept)` is for grand mean which is probably not a parameter of interest.

**4. Q: Help, what does the error message"'rank' must be a value from 'taxonomyRanks()'" mean?**

A: To start, it is recommended that you use the `taxonomyRanks(se)` function to change the taxonomy ranks in your `phyloseq` object. This is important because the rank names in your `tax_table` must be named correctly, such as "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", or "Species". If your rank names are currently labeled as something else, such as "ta1", "ta2", "ta3", and so on, you will need to update them accordingly. This issue often arises when a `tax_table` is formed from a `data.frame` instead of a `matrix`. Once you have updated the rank names, you can use the `ancombc2` function with the `tax_level` argument set to the appropriate rank level (e.g., "Genus"). This will allow you to perform statistical analyses on your microbiome data at the desired taxonomic level.

**5. Q: I ran into a problem using `rand_formula`**

A: `ancombc2` follows the `lmerTest` package in formulating the random effects. Please pay attention to the **parenthesis** and **vertical bars**. For instance, for a random subject effect, `rand_formula = "(1|subjid)" `is the correct way of specifying it, while both `rand_formula = "1|subjid"` and `rand_formula = "(subjid)"` are incorrect.

**6. Q: What are the differences between the primary results and the results of Dunnet's type of test?**

A: Say you have a `group` varible with 3 levels, `A`, `B`, and `C`, both the primary results and the results of Dunnett's type of test will provide you with differentially abundant taxa for the comparisons of `B - A` and `C - A`. However, there is a difference in the correction of p-values. The primary results only correct p-values across taxa, while Dunnett's type of test corrects p-values across taxa and for multiple comparisons (`B - A` and `C - A`), resulting in a more conservative outcome that is less prone to false positives.

**7. Q: Can `ancombc` or `ancombc2` function deal with intereaction terms?**

A: Regrettably, the inclusion of interaction terms in the `fix_formula` argument can lead to confusion in the multi-group comparisons of `ancombc2`. In this case, I suggest manually creating an interaction term to achieve the desired analysis.

**8. Q: Can you give me a more complicated example of performing ANCOM-BC2 trend test?**

A: For example, when using the trend test with a `group` variable of 5 ordered categories (`A, B, C, D, E`) in R, we are actually estimating 4 contrasts, which are (`B-A, C-A, D-A, E-A`). Testing the trend of `A < B < C < D < E` is equivalent to testing `0 < B - A < C - A < D - A < E - A`. Therefore, we can specify the contrast matrix as follows:

```R
# B-A    C-A     D-A    E-A
    1      0      0     0
    -1     1      0     0
    0     -1     1      0
    0     0     -1      1
```

In R, it should be

```R
matrix(c(1, 0, 0, 0,
  -1, 1, 0, 0,
   0, -1, 1, 0,
   0, 0, -1, 1),
   nrow = 4, 
   byrow = TRUE)
```

**9. Q: OMG, I am still very confused at structural zeros. What are they? What do `struc_zero` and `neg_lb` arguments do?**

A: A taxon is considered to have structural zeros in some (>=1) groups if it is completely or nearly completely absent in those groups. For example, if there are three groups, g1, g2, and g3, and the counts of taxon A are 0 in g1 but non-zero in g2 and g3, taxon A will be considered to contain structural zeros in g1. In this scenario, taxon A is declared to be differentially abundant between g1 and g2, g1 and g3, and is consequently globally differentially abundant with respect to the group variable. Such taxa are not further analyzed using ANCOM-BC or ANCOM-BC2, but the results are summarized in the `zero_ind`. You can treat the detection of structural zeros as performing a presence/absence test.

The detection of structural zeros is based on a separate paper, [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/). Specifically, setting `neg_lb = TRUE` indicates that both criteria stated in section 3.2 of [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/) are used to detect structural zeros. Alternatively, setting `neg_lb = FALSE` will only use equation 1 in section 3.2 of [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/) to declare structural zeros, making it a more conservative approach. As the OTU/ASV table is usually very sparse, it is recommended to choose `neg_lb = FALSE` to prevent false discoveries. However, if you have a more dense table such as a family level table with a sufficiently large sample size, using `neg_lb = TRUE` may be a better idea. It is important to note that `neg_lb` has no function if `struc_zero` is set to `FALSE`. Therefore, there are three possible combinations: `struc_zero = FALSE` (regardless of `neg_lb`), `struc_zero = TRUE, neg_lb = FALSE`, or `struc_zero = TRUE, neg_lb = TRUE`.

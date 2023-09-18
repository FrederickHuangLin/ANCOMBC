# ANCOMBC

*Author & Maintainer: Huang Lin: <huanglinfrederick@gmail.com>*

**Cautionary Notice: Please exercise discretion while using ANCOM-BC2, as it is founded on a paper that is yet to be published. The manuscript underpinning ANCOM-BC2 is presently under peer review, and as such, the existing software version may exhibit instability. The expected timeline for ANCOM-BC2's publication is set for later this year.**

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

**1. Q: What are the differences between the `formula` and `group` arguments in `ancombc` and `ancombc2`?**

A: The `formula` and `group` arguments serve different purposes in the `ancombc` and `ancombc2` functions. Here's a breakdown of their differences:

1. `formula`: This argument is used to specify the variables in your experiment that can potentially influence microbial abundances. It is essential to include all relevant variables in the `formula` to ensure proper adjustment and accurate results. For example, if you have a continuous variable like `age` as your main variable of interest, and you have additional categorical variables that need adjustment but are not directly related to your research question, you can include them in the `formula` while leaving `group` as NULL.

2. `group`: The `group` argument is optional and should only be specified if you want to detect structural zeros (presence/absence test) or perform multi-group comparisons, such as the global test, pairwise directional test, Dunnett's type of test, or trend test. If your variable of interest is a categorical variable with more than three levels and you want to conduct multi-group comparisons, you should include the `group` argument. It is important to note that `group` is not the same as `main_var` in `ancom`. In `ancombc` and `ancombc2`, `group` is used for multi-group comparisons and correction of p-values for multiple comparisons.

Remember not to include the `main_var` in the `adj_formula` in `ancom`, but always include `group` in the `formula` or `fix_formula` (in `ancombc` and `ancombc2`, respectively) if `group` is not NULL. This ensures that the appropriate adjustments and comparisons are made in the analysis.

**2. Q: Why are some taxa absent from the primary results?**

A: There are a couple of reasons why certain taxa may be absent from the primary results in `ancombc` or `ancombc2`. Here's an explanation:

1. Prevalence Exclusion: Taxa with prevalences below the specified threshold (`prv_cut`) will be excluded from the analysis. The `prv_cut` value determines the minimum prevalence required for a taxon to be considered in the analysis. If a taxon's prevalence falls below this threshold, it will not be included in the primary results.

2. Structural Zeros: Taxa that exhibit structural zeros, meaning they consistently have zero counts across all samples, will be considered significant only by the presence/absence test. The ANCOM-BC and ANCOM-BC2 methodologies are not designed to detect significant differences in taxa with structural zeros. As a result, these taxa are summarized separately and not included in the primary results of `ancombc` or `ancombc2`.

To access the results of the presence/absence test, you can refer to the `zero_ind` output. This will provide information on the taxa that exhibit structural zeros.

**3. Q: In the primary results, what do `lfc_(Intercept)`, `lfc_groupB`, and `lfc_groupC` represent if I have a `group` variable with categories `A`, `B`, and `C`?**

A: In the primary results, the terms `lfc_groupB` and `lfc_groupC` represent the log fold changes (logFC) relative to the reference group, which is group `A` by default. These logFC values indicate the difference in abundance between group `B` and group `A`, and between group `C` and group `A`, respectively. 

On the other hand, `lfc_(Intercept)` refers to the log fold change of the grand mean, which may not be a parameter of particular interest in this context.

It's worth mentioning that if you wish to change the reference group, you can use the `factor` function in R to rearrange the levels of the `group` variable accordingly.

**4. Q: I encountered an error message stating "'rank' must be a value from 'taxonomyRanks()'. What does it mean and how can I resolve it?**

A: The error message "'rank' must be a value from 'taxonomyRanks()'" typically occurs when the rank names in your `tax_table` are not properly labeled. In order to resolve this issue, it is recommended to use the `taxonomyRanks(se)` function, where `se` is your `tse` object.

Firstly, ensure that the rank names in your `tax_table` are correctly named as one of the standard taxonomic ranks, such as "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", or "Species". If the rank names are currently labeled as something else, such as "ta1", "ta2", "ta3", and so on, you will need to update them accordingly.

This issue commonly occurs when a `tax_table` is formed from a `data.frame` instead of a `matrix`. Therefore, it's important to ensure that the rank names are correctly assigned before proceeding with the analysis.

Once you have updated the rank names in your `tax_table`, you can utilize the `ancombc2` function, specifying the desired taxonomic level using the `tax_level` argument (e.g., "Genus"). This will enable you to perform statistical analyses on your microbiome data at the specified taxonomic level.

**5. Q: I encountered an issue while using `rand_formula` in `ancombc2`. What is the correct syntax for specifying random effects?**

A: When specifying random effects using `rand_formula` in `ancombc2`, it is important to follow the syntax conventions used in the `lmerTest` package. Pay close attention to the placement of **parentheses** and **vertical bars**.

To correctly specify a random subject effect, the syntax should be in the form of `"(1|subjid)"`, where `subjid` represents the variable name for the subject identifier. This syntax ensures that the random subject effect is properly accounted for in the analysis.

On the other hand, it is incorrect to use `rand_formula` as `"1|subjid"` or `"(subjid)"` for specifying random effects. The correct syntax should always include parentheses around the random effect and a vertical bar to separate it from the fixed effects.

By using the correct syntax for specifying random effects, you will be able to accurately incorporate these effects into your `ancombc2` analysis.

**6. Q: What are the differences between the primary results and the results of Dunnett's type of test in ANCOM-BC2?**

A: The primary results and the results of Dunnett's type of test in ANCOM-BC2 provide information on differentially abundant taxa, but there are differences in the correction of p-values.

In the primary results, the p-values are corrected across taxa, meaning that they account for multiple comparisons among different taxa. This correction helps control the false positive rate when determining the significance of individual taxa.

On the other hand, Dunnett's type of test not only corrects the p-values across taxa but also corrects for multiple comparisons between groups. Specifically, it compares the abundance of each taxon in groups `B` and `C` with the reference group `A`. The correction for multiple comparisons in Dunnett's type of test results in a more conservative outcome, reducing the likelihood of false positive results.

Therefore, while both the primary results and Dunnett's type of test provide information on differentially abundant taxa, the results of Dunnett's type of test offer additional control for multiple comparisons, making them more conservative and reliable.

**7. Q: Can the `ancombc` or `ancombc2` function handle interaction terms in the analysis?**

A: Unfortunately, the inclusion of interaction terms in the `fix_formula` argument of `ancombc` or `ancombc2` can lead to complexities and potential confusion in the multi-group comparisons. To address this, it is recommended to manually create the interaction term of interest outside of the formula and perform the analysis accordingly.

By manually creating the interaction term, you can ensure that the analysis accurately captures the interaction effect between variables. Once the interaction term is created, you can include it in the `fix_formula` argument or any other relevant part of the analysis, depending on your specific research question and design.

**8. Q: Can the ANCOM-BC methodology be applied to other data types such as functional abundances, RNA-seq, or single-cell RNA data?**

A: The ANCOM-BC methodology can be applied to other data types as long as they are considered compositional. However, it is essential to be aware that the methodology has been primarily benchmarked and validated using microbiome data. For more discussions, you can refer to this [post](https://github.com/FrederickHuangLin/ANCOMBC/issues/196). 

**9. Q: What does "not a positive definite matrix" mean in fitting the `ancombc2` mixed effects model? How can I debug this issue?**

A: The error message "not a positive definite matrix" indicates that the correlation matrix used in the mixed effects model is not positive definite. A positive definite matrix is a square matrix where all eigenvalues are positive. This error typically occurs when there is an issue with the data or model specification.

To debug this issue, I recommend fitting a mixed effects model to your **RAW** data using the `lmerTest` package in R. Use the same fixed effects and random effects specifications that you used in the `ancombc2` function. By fitting the model directly, you may receive more informative error messages that can help diagnose the problem.

Here are the steps you can follow to debug the issue:

1. Install and load the `lmerTest` package in R: `install.packages("lmerTest")` and `library(lmerTest)`.
2. Prepare your data in its raw format without any transformation or preprocessing.
3. Specify the fixed effects and random effects in the model formula, similar to what you used in `ancombc2`.
4. Fit the mixed effects model using the `lmer()` function from the `lmerTest` package.
5. Check if the model fitting process encounters any errors or warnings. These messages can provide valuable insights into the issue.
6. Analyze the error or warning messages to identify the underlying problem. It could be related to the data structure, model specification, or potential collinearity among variables.
7. Address the issue based on the information provided in the error or warning messages. This may involve revising the model specification, examining the data for anomalies, or resolving any collinearity issues.

By following these steps, you can gain a better understanding of the problem causing the "not a positive definite matrix" error and take appropriate actions to address it.

If you continue to encounter difficulties or need further assistance, it may be helpful to seek advice from statisticians or experts in your specific field of research.

**10. Q: If a higher LFC value corresponds to a larger abundance, why did some of my OTU/ASV counts show opposite directions?**

A: It's important to note that the log-fold change (LFC) values in the context of ANCOM-BC or ANCOM-BC2 do not directly reflect the relative abundances (such as proportions) or observed abundances (such as OTU or ASV counts). The LFC values represent the difference in bias-corrected abundances between groups.

In ANCOM-BC or ANCOM-BC2, a higher LFC value indicates a larger difference in bias-corrected abundances between groups. However, this does not necessarily mean that the group with higher LFC has a higher relative abundance or larger observed counts for a specific OTU or ASV.

**11. Q: Can you give me a more complicated example of performing ANCOM-BC2 trend test?**

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

For more in-depth discussions, you can refer to this [post](https://github.com/FrederickHuangLin/ANCOMBC/issues/204).

**12. Q: OMG, I am still very confused at structural zeros. What are they? What do the `struc_zero` and `neg_lb` arguments do?**

A: A taxon is considered to have structural zeros in some (>=1) groups if it is completely or nearly completely absent in those groups. For example, if there are three groups, g1, g2, and g3, and the counts of taxon A are 0 in g1 but non-zero in g2 and g3, taxon A will be considered to contain structural zeros in g1. In this scenario, taxon A is declared to be differentially abundant between g1 and g2, g1 and g3, and is consequently globally differentially abundant with respect to the group variable. Such taxa are not further analyzed using ANCOM-BC or ANCOM-BC2, but the results are summarized in the `zero_ind`. You can treat the detection of structural zeros as performing a presence/absence test.

The detection of structural zeros is based on a separate paper, [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/). Specifically, setting `neg_lb = TRUE` indicates that both criteria stated in section 3.2 of [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/) are used to detect structural zeros. Alternatively, setting `neg_lb = FALSE` will only use equation 1 in section 3.2 of [ANCOM-II](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/) to declare structural zeros, making it a more conservative approach. As the OTU/ASV table is usually very sparse, it is recommended to choose `neg_lb = FALSE` to prevent false discoveries. However, if you have a more dense table such as a family level table with a sufficiently large sample size, using `neg_lb = TRUE` may be a better idea. It is important to note that `neg_lb` has no function if `struc_zero` is set to `FALSE`. Therefore, there are three possible combinations: `struc_zero = FALSE` (regardless of `neg_lb`), `struc_zero = TRUE, neg_lb = FALSE`, or `struc_zero = TRUE, neg_lb = TRUE`.

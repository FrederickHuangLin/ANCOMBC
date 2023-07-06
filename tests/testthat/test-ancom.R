context("Testing ancom function")
library(ANCOMBC)
library(testthat)

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# test
test_that("`ancom` function provides expected results", {
    set.seed(123)
    out = ancom(data = tse, assay_name = "counts",
                tax_level = "Family", phyloseq = NULL,
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                main_var = "bmi_group", adj_formula = "age + nationality",
                rand_formula = NULL, lme_control = NULL,
                struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1)
    res = out$res
    test_output = res$W[1]
    expect_equal(test_output, 18)
})



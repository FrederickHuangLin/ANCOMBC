context("Testing ancombc2 function")
library(ANCOMBC)
library(testthat)

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# test
test_that("`ancombc2` function provides expected results", {
    set.seed(123)
    out = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                   fix_formula = "age + nationality + bmi_group",
                   rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = FALSE,
                   prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                   group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = 1, verbose = FALSE,
                   global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 100, verbose = FALSE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   mdfdr_control = NULL,
                   trend_control = NULL)
    res_prim = out$res
    test_output = round(res_prim$W_age[1], 2)
    expect_equal(test_output, -5.99)
})



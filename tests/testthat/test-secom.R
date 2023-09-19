context("Testing secom functions")
library(ANCOMBC)
library(testthat)

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# test 1
test_that("`secom_linear` function provides expected results", {
    set.seed(123)
    res_linear = secom_linear(data = list(tse), assay_name = "counts",
                              tax_level = "Phylum", pseudo = 0,
                              prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                              wins_quant = c(0.05, 0.95), method = "pearson",
                              soft = FALSE, thresh_len = 20, n_cv = 10,
                              thresh_hard = 0.3, max_p = 0.005, n_cl = 1)
    test_output = round(c(res_linear$corr_th[2, 1],
                          res_linear$corr_p[2, 1]), 2)
    expect_equal(test_output, c(-0.45, 0.00))
})

# test 2
test_that("`secom_dist` function provides expected results", {
    set.seed(123)
    res_dist = secom_dist(data = list(tse), assay_name = "counts",
                          tax_level = "Phylum", pseudo = 0,
                          prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                          wins_quant = c(0.05, 0.95), R = 1000,
                          thresh_hard = 0.3, max_p = 0.005, n_cl = 1)
    test_output = round(c(res_dist$dcorr_fl[2, 1],
                          res_dist$dcorr_p[2, 1]), 2)
    expect_equal(test_output, c(0.46, 0.00))
})



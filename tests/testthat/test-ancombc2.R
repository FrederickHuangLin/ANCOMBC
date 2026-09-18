context("Testing ancombc2 function")
library(ANCOMBC)
library(testthat)

data(atlas1006, package = "microbiome")

# subset to baseline
pseq = phyloseq::subset_samples(atlas1006, time == 0)

# test
test_that("`ancombc2` function provides expected results", {
    set.seed(123)
    out = ancombc2(data = pseq,
                   tax_level = "Family",
                   fix_formula = "age + nationality + bmi_group",
                   rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = FALSE,
                   prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                   group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = 1, verbose = TRUE,
                   global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 100, verbose = FALSE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   mdfdr_control = NULL,
                   trend_control = NULL)
    res_prim = out$res
    test_output = round(res_prim$W_age[1], 2)
    expect_equal(test_output, -5.99)
})

test_that("`ancombc2` function provides expected results with `conservative = FALSE`", {
    set.seed(123)
    out = ancombc2(data = pseq,
                   tax_level = "Family",
                   fix_formula = "age + nationality + bmi_group",
                   rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   conservative = FALSE,
                   prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                   group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = 1, verbose = FALSE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 100, verbose = FALSE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = NULL)

    ss_cols = c("passed_ss_age", "diff_robust_age")
    expect_true(all(ss_cols %in% colnames(out$res)))
    expect_true(all(vapply(out$res[, ss_cols], is.logical, logical(1))))

    expect_true(all(c("passed_ss", "diff_robust_abn") %in% colnames(out$res_global)))
    expect_true(all(vapply(out$res_global[, c("passed_ss", "diff_robust_abn")],
                           is.logical, logical(1))))

    for (res_tab in list(out$res_pair, out$res_dunn)) {
        ss_cols = grep("^passed_ss_|^diff_robust_", colnames(res_tab), value = TRUE)
        expect_gt(length(ss_cols), 0)
        expect_true(all(vapply(res_tab[, ss_cols], is.logical, logical(1))))
    }
})



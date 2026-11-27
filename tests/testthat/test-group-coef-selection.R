context("Testing group coefficient selection with overlapping covariate names")
library(ANCOMBC)
library(testthat)

test_that("ancombc2 multi-group tests only select coefficients of the group variable", {
    set.seed(42)
    n_tax = 60
    n_per = 12
    n_samp = 3 * n_per
    meta = data.frame(
        site = factor(rep(c("A", "B", "C"), each = n_per)),
        site_score = stats::rnorm(n_samp),
        row.names = paste0("s", seq_len(n_samp))
    )
    mu = stats::rnorm(n_tax, 3, 1)
    log_mean = mu + matrix(stats::rnorm(n_tax * n_samp, 0, 0.5), n_tax, n_samp)
    counts = matrix(stats::rpois(n_tax * n_samp, exp(log_mean)), n_tax, n_samp,
                    dimnames = list(paste0("T", seq_len(n_tax)), rownames(meta)))

    res = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta,
                 fix_formula = "site + site_score", group = "site",
                 pseudo_sens = FALSE, prv_cut = 0,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 20),
                 trend_control = list(
                     contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)),
                     node = list(2),
                     B = 20
                 ),
                 verbose = FALSE)
    ))

    # Pairwise and Dunnett results must only contain siteB and siteC, not site_score
    expect_false(any(grepl("site_score", colnames(res$res_pair), fixed = TRUE)))
    expect_false(any(grepl("site_score", colnames(res$res_dunn), fixed = TRUE)))
    expect_identical(
        grep("^lfc_", colnames(res$res_pair), value = TRUE),
        c("lfc_siteB", "lfc_siteC", "lfc_siteC_siteB")
    )
    expect_identical(
        grep("^lfc_", colnames(res$res_dunn), value = TRUE),
        c("lfc_siteB", "lfc_siteC")
    )
    expect_identical(
        grep("^lfc_", colnames(res$res_trend), value = TRUE),
        c("lfc_siteB", "lfc_siteC")
    )

    # Selected pairs must also work when another covariate starts with group's name
    res_sel = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta,
                 fix_formula = "site + site_score", group = "site",
                 pseudo_sens = FALSE, prv_cut = 0,
                 pairwise = TRUE,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 20,
                                      pairs = list(c("B", "C"))),
                 verbose = FALSE)
    ))
    expect_identical(grep("^lfc_", colnames(res_sel$res_pair), value = TRUE), "lfc_siteC_siteB")

    # Global test must match running with a renamed non-overlapping covariate
    meta_renamed = meta
    colnames(meta_renamed)[2] = "cov_score"
    res_ref = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta_renamed,
                 fix_formula = "site + cov_score", group = "site",
                 pseudo_sens = FALSE, prv_cut = 0,
                 global = TRUE, verbose = FALSE)
    ))
    expect_equal(res$res_global$W, res_ref$res_global$W, tolerance = 1e-10)
    expect_equal(res$res_global$p_val, res_ref$res_global$p_val, tolerance = 1e-10)
})

test_that("ancombc structural zero masking and global test ignore covariates containing group name", {
    set.seed(43)
    n_tax = 60
    n_per = 12
    n_samp = 3 * n_per
    meta = data.frame(
        site = factor(rep(c("A", "B", "C"), each = n_per)),
        site_score = stats::rnorm(n_samp),
        row.names = paste0("s", seq_len(n_samp))
    )
    counts = matrix(stats::rpois(n_tax * n_samp, lambda = 20), n_tax, n_samp,
                    dimnames = list(paste0("T", seq_len(n_tax)), rownames(meta)))
    counts[1, meta$site == "B"] = 0

    meta_renamed = meta
    colnames(meta_renamed)[2] = "cov_score"

    res = suppressWarnings(suppressMessages(
        ancombc(data = counts, meta_data = meta,
                formula = "site + site_score", group = "site",
                struc_zero = TRUE, global = TRUE, prv_cut = 0, verbose = FALSE)
    ))
    res_ref = suppressWarnings(suppressMessages(
        ancombc(data = counts, meta_data = meta_renamed,
                formula = "site + cov_score", group = "site",
                struc_zero = TRUE, global = TRUE, prv_cut = 0, verbose = FALSE)
    ))

    expect_equal(res$res$se$site_score, res_ref$res$se$cov_score, tolerance = 1e-10)
    expect_equal(res$res$p_val$site_score, res_ref$res$p_val$cov_score, tolerance = 1e-10)
    expect_equal(res$res_global$W, res_ref$res_global$W, tolerance = 1e-10)

    # Non-conservative sensitivity analysis with pairwise comparisons also ignores site_score
    res_sens = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta,
                 fix_formula = "site + site_score", group = "site",
                 pseudo_sens = TRUE, conservative = FALSE,
                 prv_cut = 0, pairwise = TRUE, global = TRUE,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 20),
                 verbose = FALSE)
    ))
    expect_false(any(grepl("site_score", colnames(res_sens$res_pair), fixed = TRUE)))
})

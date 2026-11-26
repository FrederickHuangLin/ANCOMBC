context("Testing estimators, tests, and options of method_control")
library(ANCOMBC)
library(testthat)

# Poisson log-normal counts: n_tax taxa, three groups of n_per samples, a
# continuous covariate, and taxa 1 to n_da with a group B effect
sim_counts = function(n_tax = 60, n_per = 15, n_da = 6, seed = 1) {
    set.seed(seed)
    n = 3 * n_per
    meta = data.frame(group = factor(rep(c("A", "B", "C"), each = n_per)),
                      age = round(stats::rnorm(n, 50, 10)),
                      row.names = paste0("s", seq_len(n)))
    mu = stats::rnorm(n_tax, 3, 1)
    lfc = matrix(0, n_tax, 3)
    lfc[seq_len(n_da), 2] = 1.5
    log_mean = mu + lfc[, as.integer(meta$group)] +
        matrix(stats::rnorm(n_tax * n, 0, 0.5), n_tax, n)
    counts = matrix(stats::rpois(n_tax * n, exp(log_mean)), n_tax, n,
                    dimnames = list(paste0("T", seq_len(n_tax)),
                                    rownames(meta)))
    list(counts = counts, meta = meta)
}

test_that("the sandwich variance does not depend on the location of a covariate", {
    set.seed(1)
    n = 150
    grp = rep(0:1, length.out = n)
    age = stats::rnorm(n, 5000, 1000)
    eps = matrix(stats::rnorm(5 * n), 5, n)
    x1 = cbind(1, grp, age)
    x2 = cbind(1, grp, age - 5000)
    v1 = ANCOMBC:::.sandwich_vcov(x1, eps, c("int", "grp", "age"))$var_hat
    v2 = ANCOMBC:::.sandwich_vcov(x2, eps, c("int", "grp", "age"))$var_hat
    expect_equal(v1[, 2], v2[, 2], tolerance = 1e-8)
    expect_equal(v1[, 3], v2[, 3], tolerance = 1e-8)
})

test_that("vcov = 'observed' is the HC0 estimator of the samples of each taxon", {
    set.seed(2)
    n = 40
    x = cbind(1, rep(0:1, each = n / 2), stats::rnorm(n))
    eps = matrix(stats::rnorm(3 * n), 3, n)
    eps[1, c(1, 5, 30)] = NA
    eps[3, x[, 2] == 0] = NA
    fix_eff = c("int", "grp", "z")
    out = ANCOMBC:::.sandwich_vcov(x, eps, fix_eff, vcov = "observed")
    rows = !is.na(eps[1, ])
    xr = x[rows, ]
    bread = solve(crossprod(xr))
    hc0 = bread %*% crossprod(xr * eps[1, rows]) %*% bread
    expect_equal(unname(out$vcov_hat[[1]]), unname(hc0), tolerance = 1e-10)
    # Taxon 3 is absent from the reference level: the intercept and the group
    # coefficient are not estimable
    expect_true(all(is.na(out$var_hat[3, 1:2])))
    expect_false(is.na(out$var_hat[3, 3]))
})

test_that("non-estimable coefficients are identified from the samples of each taxon", {
    x = cbind(1, rep(c(0, 1, 0), each = 4), rep(c(0, 0, 1), each = 4))
    colnames(x) = c("(Intercept)", "grpB", "grpC")
    y = matrix(stats::rnorm(24), 2, 12)
    y[2, 1:4] = NA
    est = ANCOMBC:::.estimable_mask(x, y)
    expect_true(all(est[1, ]))
    expect_equal(unname(est[2, ]), c(FALSE, FALSE, FALSE))
})

test_that("ancombc2 reports NA for a taxon absent from the reference level", {
    sim = sim_counts(seed = 3)
    sim$counts[60, sim$meta$group == "A"] = 0
    out = suppressWarnings(suppressMessages(
        ancombc2(data = sim$counts, meta_data = sim$meta,
                 fix_formula = "group", group = "group", pseudo_sens = FALSE,
                 prv_cut = 0, verbose = FALSE)))
    res = out$res[out$res$taxon == "T60", ]
    expect_true(is.na(res$lfc_groupB))
    expect_true(is.na(res$lfc_groupC))
    expect_false(res$diff_groupB)
    expect_false(res$diff_groupC)
})

test_that("the variance update of the E-M algorithm attains an interior minimum", {
    kappa = ANCOMBC:::.bias_em_kappa(beta = 1, nu0 = 0.99, weights = 1,
                                     centre = 0, current = 1)
    expect_equal(kappa, 0.01, tolerance = 1e-6)
    kappa0 = ANCOMBC:::.bias_em_kappa(beta = 0.5, nu0 = 1, weights = 1,
                                      centre = 0, current = 1)
    expect_equal(kappa0, 0)
})

test_that("global_p = 'upper' gives the upper tail p-value of the Wald test", {
    x = cbind(1, diag(3)[rep(1:3, each = 2), 2:3])
    colnames(x) = c("(Intercept)", "groupB", "groupC")
    beta_hat = matrix(c(0.1, 0.1, 2, 1), 2, 2, byrow = TRUE,
                      dimnames = list(c("t1", "t2"), c("groupB", "groupC")))
    beta_hat = cbind("(Intercept)" = 0, beta_hat)
    vcov_hat = rep(list(diag(3)), 2)
    dof = matrix(20, 2, 3)
    up = ANCOMBC:::.ancombc_global_F(x, "group", beta_hat, vcov_hat, dof = dof,
                                     p_adj_method = "none", alpha = 0.05,
                                     global_p = "upper")
    W = c(0.02, 5)
    expect_equal(up$W, W)
    expect_equal(up$p_val, stats::pf(W / 2, 2, 20, lower.tail = FALSE))
    two = ANCOMBC:::.ancombc_global_F(x, "group", beta_hat, vcov_hat, dof = dof,
                                      p_adj_method = "none", alpha = 0.05)
    expect_equal(two$p_val, 2 * pmin(stats::pf(W, 2, 20),
                                      stats::pf(W, 2, 20, lower.tail = FALSE)))
})

test_that("dunnet_null = 'estimable' changes only taxa with a missing statistic", {
    x = cbind(1, diag(4)[rep(1:4, each = 2), 2:4])
    colnames(x) = c("(Intercept)", "groupB", "groupC", "groupD")
    W = matrix(c(3, 0.5, 0.2, 3, NA, NA), 2, 3, byrow = TRUE,
               dimnames = list(c("t1", "t2"), colnames(x)[2:4]))
    dof = matrix(30, 2, 3)
    set.seed(4)
    all_cells = ANCOMBC:::.dunn_global(x, "group", W, B = 500, dof = dof,
                                       p_adj_method = "none", alpha = 0.05)
    set.seed(4)
    est_cells = ANCOMBC:::.dunn_global(x, "group", W, B = 500, dof = dof,
                                       p_adj_method = "none", alpha = 0.05,
                                       dunnet_null = "estimable")
    expect_equal(all_cells$p_val[1], est_cells$p_val[1])
    expect_lt(est_cells$p_val[2], all_cells$p_val[2])
})

test_that("ancom tests a multi-level main variable given adj_formula", {
    sim = sim_counts(n_tax = 12, seed = 5)
    out = suppressWarnings(suppressMessages(
        ancom(data = sim$counts, meta_data = sim$meta, main_var = "group",
              adj_formula = "age", prv_cut = 0, verbose = FALSE)))
    comp = log(sim$counts + 1)
    df = data.frame(x = comp[2, ] - comp[1, ], sim$meta)
    p_partial = stats::anova(stats::lm(x ~ age + group, data = df))["group", "Pr(>F)"]
    expect_equal(out$p_data["T1", "T2"], p_partial, tolerance = 1e-10)
})

test_that("the global test of the sensitivity analysis is adjusted for other terms", {
    set.seed(6)
    meta = data.frame(group = factor(rep(c("A", "B", "C"), each = 10)),
                      age = stats::rnorm(30))
    y = stats::rnorm(30) + meta$age
    p = ANCOMBC:::.ancombc2_sens_p(y, meta_data = meta,
                                   fix_formula = "group + age",
                                   fix_eff = c("(Intercept)", "groupB",
                                               "groupC", "age"),
                                   group = "group", group_eff = NULL,
                                   pairwise = FALSE, global = TRUE)
    p_partial = stats::anova(stats::lm(y ~ age + group, data = meta))["group", "Pr(>F)"]
    expect_equal(unname(p["global"]), p_partial, tolerance = 1e-10)
})

test_that("method_control is validated", {
    sim = sim_counts(seed = 7)
    expect_error(ancombc2(data = sim$counts, meta_data = sim$meta,
                          fix_formula = "group", group = "group",
                          method_control = list(vcov = "hc3"), verbose = FALSE),
                 "method_control\\$vcov")
    expect_error(ancombc2(data = sim$counts, meta_data = sim$meta,
                          fix_formula = "group", group = "group",
                          method_control = list(variance = "observed"),
                          verbose = FALSE),
                 "Unknown entries")
})

test_that("the options of method_control run and keep the output structure", {
    sim = sim_counts(seed = 8)
    run = function(ctrl) suppressWarnings(suppressMessages(
        ancombc2(data = sim$counts, meta_data = sim$meta,
                 fix_formula = "group + age", group = "group",
                 pseudo_sens = FALSE, prv_cut = 0, global = TRUE,
                 dunnet = TRUE, verbose = FALSE,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 20),
                 method_control = ctrl)))
    set.seed(9)
    base = run(list())
    set.seed(9)
    alt = run(list(vcov = "observed", global_p = "upper",
                   failed_taxa = "exclude", dunnet_null = "estimable"))
    expect_identical(colnames(base$res), colnames(alt$res))
    expect_identical(colnames(base$res_global), colnames(alt$res_global))
    expect_identical(colnames(base$res_dunn), colnames(alt$res_dunn))
    # No taxon fails entirely, so the sampling fractions agree up to the
    # variance estimator used by the E-M algorithm
    expect_equal(length(base$samp_frac), length(alt$samp_frac))
    expect_true(all(alt$res_global$p_val >= 0 & alt$res_global$p_val <= 1))
})

test_that("pairwise comparisons can be restricted to selected pairs", {
    sim = sim_counts(seed = 10)
    run = function(pairs) suppressWarnings(suppressMessages(
        ancombc2(data = sim$counts, meta_data = sim$meta,
                 fix_formula = "group", group = "group", pseudo_sens = FALSE,
                 prv_cut = 0, pairwise = TRUE, verbose = FALSE,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 20,
                                      pairs = pairs))))
    full = run(NULL)
    sel = run(list(c("C", "B"), c("A", "B")))
    expect_identical(colnames(sel$res_pair),
                     c("taxon", "lfc_groupC_groupB", "lfc_groupB",
                       "se_groupC_groupB", "se_groupB",
                       "W_groupC_groupB", "W_groupB",
                       "p_groupC_groupB", "p_groupB",
                       "q_groupC_groupB", "q_groupB",
                       "diff_groupC_groupB", "diff_groupB"))
    expect_equal(sel$res_pair$lfc_groupC_groupB, full$res_pair$lfc_groupC_groupB)
    expect_equal(sel$res_pair$p_groupB, full$res_pair$p_groupB)
    # Fewer comparisons are counted in the mdFDR adjustment
    expect_true(all(sel$res_pair$q_groupB <= full$res_pair$q_groupB + 1e-12))
    expect_error(run(list(c("A", "D"))), "mdfdr_control\\$pairs")
})

test_that("variable names that are not syntactically valid are supported", {
    sim = sim_counts(seed = 11)
    meta = sim$meta
    colnames(meta) = c("body-site", "age (years)")
    out = suppressWarnings(suppressMessages(
        ancombc2(data = sim$counts, meta_data = meta,
                 fix_formula = "`body-site` + `age (years)`",
                 group = "body-site", pseudo_sens = FALSE, prv_cut = 0,
                 global = TRUE, verbose = FALSE)))
    ref = suppressWarnings(suppressMessages(
        ancombc2(data = sim$counts, meta_data = sim$meta,
                 fix_formula = "group + age", group = "group",
                 pseudo_sens = FALSE, prv_cut = 0, global = TRUE,
                 verbose = FALSE)))
    expect_true("lfc_`body-site`B" %in% colnames(out$res))
    expect_equal(out$res[["lfc_`body-site`B"]], ref$res$lfc_groupB)
    expect_equal(out$res_global$p_val, ref$res_global$p_val)

    out_bc = suppressWarnings(suppressMessages(
        ancombc(data = sim$counts, meta_data = meta,
                formula = "`body-site` + `age (years)`", group = "body-site",
                prv_cut = 0, verbose = FALSE)))
    expect_true("`body-site`B" %in% colnames(out_bc$res$lfc))
})

context("Testing bug fixes")
library(ANCOMBC)
library(testthat)

# Negative binomial counts with additional zeros
sim_counts = function(n_tax = 40, n_samp = 24, seed = 1, zero_p = 0.1) {
    set.seed(seed)
    mu = exp(rnorm(n_tax, 4, 1))
    counts = t(vapply(mu, function(m) rnbinom(n_samp, mu = m, size = 2),
                      FUN.VALUE = numeric(n_samp)))
    counts[matrix(runif(n_tax * n_samp) < zero_p, n_tax)] = 0
    rownames(counts) = paste0("T", seq_len(n_tax))
    colnames(counts) = paste0("S", seq_len(n_samp))
    counts
}

test_that("a failed mixed-model structure for the last pattern gives NA", {
    counts = sim_counts(seed = 11)
    meta = data.frame(grp = rep(c("A", "B"), each = 12),
                      subj = rep(paste0("s", 1:6), 4),
                      row.names = colnames(counts))
    # The last taxon has a single level of grp among its nonzero samples
    counts[40, meta$grp == "B"] = 0
    out = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta, fix_formula = "grp",
                 rand_formula = "(1 | subj)", prv_cut = 0,
                 pseudo_sens = FALSE, verbose = FALSE)))
    expect_equal(nrow(out$res), 40)
    expect_true(all(is.na(out$res[40, c("lfc_(Intercept)", "lfc_grpB")])))
    expect_true(all(is.finite(out$res$lfc_grpB[1:39])))
})

test_that("a taxon with a failed fit is not rejected by the Dunnett screen", {
    x = matrix(0, nrow = 1, ncol = 4,
               dimnames = list(NULL, c("(Intercept)", "age",
                                       "batchb2", "batchb3")))
    set.seed(1)
    W = matrix(c(rnorm(18), NA, NA), nrow = 10, ncol = 2, byrow = TRUE)
    dof = matrix(20, nrow = 10, ncol = 2)
    # A failed mixed-model fit has missing degrees of freedom
    dof[10, ] = NA
    res = expect_silent(ANCOMBC:::.dunn_global(x = x, group = "batch", W = W,
                                               B = 100, dof = dof,
                                               p_adj_method = "BH",
                                               alpha = 0.05))
    expect_equal(res$p_val[10], 1)
    expect_false(res$diff_abn[10])
    # The other taxa use the same random draws as with complete df
    set.seed(1)
    W = matrix(c(rnorm(18), NA, NA), nrow = 10, ncol = 2, byrow = TRUE)
    res_ref = ANCOMBC:::.dunn_global(x = x, group = "batch", W = W[1:9, ],
                                     B = 100, dof = dof[1:9, ],
                                     p_adj_method = "BH", alpha = 0.05)
    expect_equal(res$p_val[1:9], res_ref$p_val)
})

test_that("the mixed-model check does not depend on the first taxon", {
    counts = sim_counts(n_tax = 60, n_samp = 40, seed = 2, zero_p = 0.2)
    meta = data.frame(body.site = rep(c("gut", "skin"), each = 20),
                      subject = rep(paste0("s", 1:4), times = 10),
                      row.names = colnames(counts))
    # The first taxon is observed in a single body site
    counts[1, meta$body.site == "skin"] = 0
    counts[1, meta$body.site == "gut"] = pmax(counts[1, meta$body.site == "gut"], 5)
    fit = function(cnt) {
        suppressWarnings(suppressMessages(
            ancombc2(data = cnt, meta_data = meta, fix_formula = "body.site",
                     rand_formula = "(1 | subject)", pseudo_sens = FALSE,
                     verbose = FALSE)))
    }
    out1 = fit(counts)
    out2 = fit(counts[c(2, 1, 3:60), ])
    expect_true(all(is.na(out1$res[out1$res$taxon == "T1", "lfc_body.siteskin"])))
    res1 = out1$res[match(out2$res$taxon, out1$res$taxon), ]
    rownames(res1) = NULL
    expect_equal(res1, out2$res)
})

test_that("the random-effects global test allows group as the first term", {
    counts = sim_counts(n_tax = 40, n_samp = 36, seed = 8)
    meta = data.frame(grp = rep(c("A", "B", "C"), each = 12),
                      subj = rep(paste0("s", 1:6), 6),
                      x = rnorm(36), row.names = colnames(counts))
    counts[1:6, meta$grp == "C"] = counts[1:6, meta$grp == "C"] * 10
    fit = function(f) {
        suppressWarnings(suppressMessages(
            ancombc2(data = counts, meta_data = meta, fix_formula = f,
                     rand_formula = "(1 | subj)", prv_cut = 0,
                     group = "grp", global = TRUE, pseudo_sens = FALSE,
                     verbose = FALSE)))
    }
    out1 = fit("grp")
    expect_false(anyNA(out1$res_global$p_val))
    expect_true(all(out1$res_global$W > 0))
    expect_true(any(out1$res_global$diff_abn))
    out2 = fit("grp + x")
    out3 = fit("x + grp")
    expect_equal(out2$res_global, out3$res_global)
})

test_that("messages that list names accept '%' and long lists", {
    cnt = matrix(1:60, nrow = 6, ncol = 10,
                 dimnames = list(paste0("T", 1:6), paste0("S", 1:10)))
    meta = data.frame(g = rep(c("50%s", "b%d", "c"), c(3, 3, 4)),
                      row.names = colnames(cnt))
    expect_warning(data_sanity_check(cnt, meta_data = meta, fix_formula = "g",
                                     group = "g", verbose = FALSE),
                   "50%s, b%d, c", fixed = TRUE)
    meta$g[1:2] = "c"
    expect_error(suppressWarnings(
        data_sanity_check(cnt, meta_data = meta, fix_formula = "g",
                          group = "g", verbose = FALSE)),
        "following group(s): \n50%s", fixed = TRUE)

    # Samples without nonzero counts, with long names containing '%'
    set.seed(1)
    cnt = matrix(rnbinom(100 * 400, mu = 50, size = 1), nrow = 100,
                 dimnames = list(paste0("T", 1:100),
                                 paste0("sample_%d_with_a_long_name_",
                                        sprintf("%05d", 1:400))))
    cnt[, 1:300] = 0
    meta = data.frame(grp = rep(c("a", "b"), length.out = 400),
                      row.names = colnames(cnt))
    expect_warning(suppressMessages(
        ancombc2(data = cnt, meta_data = meta, fix_formula = "grp",
                 prv_cut = 0.1, lib_cut = 0, pseudo_sens = FALSE,
                 verbose = FALSE)),
        "sample_%d_with_a_long_name_00010, ... and 290 more", fixed = TRUE)
})

test_that("the aliases rank and assay.type take effect", {
    counts = sim_counts(n_tax = 60, n_samp = 30, seed = 3)
    tax = data.frame(Kingdom = "Bacteria",
                     Phylum = paste0("P", rep(1:3, each = 20)),
                     Genus = paste0("G", rep(1:12, each = 5)),
                     row.names = rownames(counts))
    meta = data.frame(grp = rep(c("A", "B", "C"), each = 10),
                      row.names = colnames(counts))
    ps = phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
                            phyloseq::tax_table(as.matrix(tax)),
                            phyloseq::sample_data(meta))
    tse = TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(cnt = counts), rowData = S4Vectors::DataFrame(tax),
        colData = S4Vectors::DataFrame(meta))
    quiet = function(expr) suppressWarnings(suppressMessages(expr))

    out_alias = quiet(ancombc2(ps, rank = "Genus", fix_formula = "grp",
                               pseudo_sens = FALSE, verbose = FALSE))
    out_main = quiet(ancombc2(ps, tax_level = "Genus", fix_formula = "grp",
                              pseudo_sens = FALSE, verbose = FALSE))
    expect_equal(nrow(out_alias$res), 12)
    expect_equal(out_alias$res, out_main$res)
    out = quiet(ancombc2(tse, assay.type = "cnt", rank = "Genus",
                         fix_formula = "grp", pseudo_sens = FALSE,
                         verbose = FALSE))
    expect_equal(sort(out$res$taxon), sort(out_main$res$taxon))
    out = quiet(ancombc(tse, assay.type = "cnt", rank = "Genus",
                        formula = "grp", verbose = FALSE))
    expect_equal(nrow(out$feature_table), 12)
    out = quiet(ancom(ps, rank = "Genus", main_var = "grp", verbose = FALSE))
    expect_equal(nrow(out$res), 12)
    out = quiet(secom_linear(list(ps), rank = "Genus", prv_cut = 0,
                             lib_cut = 0, thresh_len = 5, n_cv = 2))
    expect_equal(nrow(out$corr_th), 12)
    out = quiet(secom_dist(list(ps), rank = "Genus", prv_cut = 0,
                           lib_cut = 0, R = 10))
    expect_equal(nrow(out$dcorr_fl), 12)
})

test_that("SECOM uses the same samples for the original and aggregated tables", {
    set.seed(4)
    asv = matrix(rnbinom(120 * 16, mu = 40, size = 1), nrow = 120, ncol = 16,
                 dimnames = list(paste0("ASV", 1:120), paste0("S", 1:16)))
    genus = paste0("G", rep(1:60, each = 2))
    # ASV1 is observed in S1 only and holds most of its reads, so S1 passes
    # lib_cut at the genus level but not after the prevalence filter of ASVs
    asv["ASV1", ] = 0
    asv["ASV1", "S1"] = 5000
    asv[-1, "S1"] = rpois(119, 5)
    agg = rowsum(asv, genus)
    meta = data.frame(x = rnorm(16), row.names = colnames(asv))
    out = suppressWarnings(suppressMessages(
        secom_linear(data = list(asv), aggregate_data = list(agg),
                     meta_data = list(meta), prv_cut = 0.5, lib_cut = 1000,
                     verbose = FALSE, n_cv = 2, thresh_len = 5)))
    expect_equal(ncol(out$y_hat), 15)
    expect_false("S1" %in% colnames(out$y_hat))
})

test_that("an informative error is given when no taxon can be fitted", {
    counts = sim_counts(n_tax = 60, n_samp = 30, seed = 4, zero_p = 0)
    meta = data.frame(batch = rep(c("B1", "B2", "B3"), each = 10),
                      row.names = colnames(counts))
    # Every taxon is observed in a single batch
    tax_batch = rep(c("B1", "B2", "B3"), 20)
    for (i in seq_len(nrow(counts))) counts[i, meta$batch != tax_batch[i]] = 0
    expect_error(suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta, fix_formula = "batch",
                 prv_cut = 0, pseudo_sens = FALSE, verbose = FALSE))),
        "cannot be estimated for the following coefficients:\n(Intercept), batchB2, batchB3",
        fixed = TRUE)
})

test_that("the cluster is stopped when a parallel run fails", {
    counts = sim_counts(n_tax = 20, n_samp = 12, seed = 5)
    meta = data.frame(grp = rep(c("A", "B"), each = 6),
                      row.names = colnames(counts))
    n_con = nrow(showConnections())
    # The formula refers to a variable that is not in the metadata
    expect_error(suppressMessages(
        ancombc2(data = counts, meta_data = meta, fix_formula = "grp + age",
                 pseudo_sens = FALSE, n_cl = 2, verbose = FALSE)),
        "not in the meta data")
    expect_equal(nrow(showConnections()), n_con)
    expect_equal(foreach::getDoParName(), "doSEQ")
    expect_error(suppressMessages(
        ancombc(data = counts, meta_data = meta, formula = "grp + age",
                n_cl = 2, verbose = FALSE)),
        "not in the meta data")
    expect_equal(nrow(showConnections()), n_con)
    expect_equal(foreach::getDoParName(), "doSEQ")
})

test_that("a clear error is given when group is not in the formula", {
    counts = sim_counts(n_tax = 30, n_samp = 30, seed = 9)
    meta = data.frame(grp = rep(c("A", "B", "C"), each = 10), x = rnorm(30),
                      row.names = colnames(counts))
    expect_error(suppressMessages(
        ancombc2(data = counts, meta_data = meta, fix_formula = "x",
                 group = "grp", dunnet = TRUE, pseudo_sens = FALSE,
                 verbose = FALSE)),
        "The group variable `grp` is not included in the formula",
        fixed = TRUE)
    expect_error(suppressMessages(
        data_sanity_check(data = counts, meta_data = meta, fix_formula = "x",
                          group = "grp", global = TRUE, verbose = FALSE)),
        "please add `grp` to the formula", fixed = TRUE)
    # Detection of structural zeros does not require group in the formula
    res = suppressMessages(
        data_sanity_check(data = counts, meta_data = meta, fix_formula = "x",
                          group = "grp", struc_zero = TRUE, verbose = FALSE))
    expect_equal(dim(res$feature_table), c(30L, 30L))
})

test_that("taxa with zero variances are excluded from the E-M estimation", {
    set.seed(12)
    beta = c(rnorm(50, 0.5, 0.2), 2)
    var_hat = c(runif(50, 0.01, 0.05), 0)
    names(beta) = names(var_hat) = paste0("T", 1:51)
    # The taxon with a zero variance does not contribute to the estimate
    expect_equal(ANCOMBC:::.bias_em(beta, var_hat, tol = 1e-5, max_iter = 100),
                 ANCOMBC:::.bias_em(beta[1:50], var_hat[1:50],
                                    tol = 1e-5, max_iter = 100))
    beta_mat = cbind(a = beta, b = beta)
    var_mat = cbind(a = var_hat, b = var_hat + 0.01)
    expect_message(ANCOMBC:::.check_bias_input(beta_mat, var_mat),
                   "following taxa:\nT51\nThese taxa are excluded from the estimation of the sample-specific bias for: a",
                   fixed = TRUE)
})

test_that("non-estimable coefficients are NA and are not tested", {
    counts = sim_counts(n_tax = 60, n_samp = 30, seed = 3)
    meta = data.frame(grp = rep(c("A", "B", "C"), each = 10),
                      subj = rep(paste0("s", 1:5), 6), x = rnorm(30),
                      row.names = colnames(counts))
    # T1 and T2 are absent from group C; T1 differs between A and B
    counts[1:2, meta$grp == "C"] = 0
    counts[1, meta$grp == "B"] = counts[1, meta$grp == "B"] * 20
    tc = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                              matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)),
              node = list(2, 2), B = 20)
    fit = function(rand_formula) {
        set.seed(1)
        suppressWarnings(suppressMessages(
            ancombc2(data = counts, meta_data = meta, fix_formula = "x + grp",
                     rand_formula = rand_formula, group = "grp", prv_cut = 0,
                     global = TRUE, pairwise = TRUE, dunnet = TRUE,
                     trend = TRUE, trend_control = tc, pseudo_sens = FALSE,
                     verbose = FALSE,
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 20))))
    }
    out = fit(NULL)
    expect_true(all(is.na(out$res[1:2, c("lfc_grpC", "se_grpC", "W_grpC")])))
    expect_equal(out$res$p_grpC[1:2], c(1, 1))
    expect_false(any(out$res$diff_grpC[1:2]))
    expect_true(all(is.finite(out$res$lfc_grpB[1:2])))
    expect_true(all(is.na(out$res_pair$lfc_grpC_grpB[1:2])))
    expect_false(any(out$res_pair$diff_grpC_grpB[1:2]))
    expect_true(all(is.na(out$res_dunn$lfc_grpC[1:2])))
    expect_false(any(out$res_dunn$diff_grpC[1:2]))
    # The global test uses the estimable group coefficient
    expect_true(out$res_global$diff_abn[1])
    expect_equal(out$res_trend$p_val[1:2], c(1, 1))
    expect_false(any(out$res_trend$diff_abn[1:2]))
    expect_true(all(is.finite(out$res$lfc_grpC[-(1:2)])))

    out = fit("(1 | subj)")
    expect_true(all(is.na(out$res$lfc_grpC[1:2])))
    expect_true(all(is.na(out$res_pair$lfc_grpC_grpB[1:2])))
    expect_equal(out$res_trend$p_val[1:2], c(1, 1))

    # A taxon absent from the reference group has a non-estimable coefficient
    # whose column is 1 in samples where the taxon is observed. Its bias
    # enters the sampling fractions as for a coefficient of 0.
    counts = sim_counts(n_tax = 60, n_samp = 30, seed = 3)
    counts[3, meta$grp == "A"] = 0
    out = suppressWarnings(suppressMessages(
        ancombc2(data = counts, meta_data = meta, fix_formula = "x + grp",
                 prv_cut = 0, pseudo_sens = FALSE, verbose = FALSE)))
    expect_true(anyNA(out$res[3, c("lfc_grpB", "lfc_grpC")]))
    o1 = log(counts)
    o1[is.infinite(o1)] = NA
    y1 = o1 - rowMeans(o1, na.rm = TRUE)
    x = stats::model.matrix(~ x + grp, data = meta)
    beta1 = ANCOMBC:::.iter_mle(x = x, y = y1, meta_data = meta,
                                formula = "x + grp", theta = NULL,
                                tol = 1e-2, max_iter = 20)$beta
    beta1[is.na(beta1)] = 0
    beta1 = t(t(beta1) - out$delta_em)
    theta = colMeans(y1 - beta1 %*% t(x), na.rm = TRUE)
    expect_equal(out$samp_frac, theta, tolerance = 1e-10)
})

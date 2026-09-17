# ---------------------------------------------------------------------------
# edge_cases.R
#
# Usage: Rscript edge_cases.R <libpath> <outfile.rds>
#
# Minimal reproductions of the edge cases that the vignette workloads do not
# reach: non-finite responses, per-taxon rank deficiency, an all-zero taxon,
# a data-dependent formula term, an unused level of the group factor, and a
# single-level group. Each case returns either its result or the error message,
# so the two libraries can be compared without either run aborting.
# ---------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
libpath = normalizePath(args[1], mustWork = TRUE)
outfile = args[2]
.libPaths(c(libpath, .libPaths()))
suppressPackageStartupMessages(library(ANCOMBC, lib.loc = libpath))
cat(sprintf("ANCOMBC %s from %s\n", as.character(packageVersion("ANCOMBC")),
            getNamespaceInfo("ANCOMBC", "path")))

# Synthetic count table: 60 taxa x 60 samples, three balanced groups.
.make_data = function(seed = 1, n = 60, d = 60) {
    set.seed(seed)
    cnt = matrix(rnbinom(d * n, mu = 200, size = 2), nrow = d, ncol = n)
    rownames(cnt) = paste0("T", seq_len(d))
    colnames(cnt) = paste0("S", seq_len(n))
    md = data.frame(cat_cov = factor(rep(c("a", "b", "c"), each = n / 3)),
                    cont_cov = rnorm(n),
                    row.names = colnames(cnt))
    list(cnt = cnt, md = md)
}

.try = function(expr) {
    r = try(suppressWarnings(suppressMessages(eval(expr))), silent = TRUE)
    if (inherits(r, "try-error"))
        list(status = "error", msg = conditionMessage(attr(r, "condition")))
    else list(status = "ok", value = r)
}

out = list()

## Case 1: ancombc() with pseudo = 0 on a table containing zeros.
## log(0) = -Inf enters the response of the fixed-effects fit.
d1 = .make_data()
d1$cnt[1, 1:5] = 0
d1$cnt[2, 3:9] = 0
out$ancombc_pseudo0 = .try(quote(
    ancombc(data = d1$cnt, meta_data = d1$md, formula = "cont_cov + cat_cov",
            pseudo = 0, p_adj_method = "holm", prv_cut = 0, lib_cut = 0,
            group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
            tol = 1e-5, max_iter = 20, conserve = TRUE, alpha = 0.05,
            global = TRUE, n_cl = 1, verbose = FALSE)))

## Case 2: a taxon whose observed samples cover only one level of cat_cov.
## Its usable-sample design is rank deficient while the shared design is not.
d2 = .make_data(seed = 2)
d2$cnt[1, d2$md$cat_cov != "a"] = 0
d2$cnt[2, d2$md$cat_cov == "c"] = 0
out$ancombc2_rank_deficient = .try(quote(
    ancombc2(data = d2$cnt, meta_data = d2$md,
             fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
             p_adj_method = "holm", pseudo_sens = FALSE,
             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
             group = "cat_cov", struc_zero = TRUE, neg_lb = TRUE,
             alpha = 0.05, n_cl = 1, verbose = FALSE,
             global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)))

## Case 3: an all-zero taxon, whose log-abundance row is entirely missing.
d3 = .make_data(seed = 3)
d3$cnt[1, ] = 0
out$ancombc2_all_zero_taxon = .try(quote(
    ancombc2(data = d3$cnt, meta_data = d3$md,
             fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
             p_adj_method = "holm", pseudo_sens = FALSE,
             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
             group = "cat_cov", struc_zero = TRUE, neg_lb = TRUE,
             alpha = 0.05, n_cl = 1, verbose = FALSE,
             global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)))

## Case 4: a data-dependent term. poly() rebuilds its basis from whichever
## samples a taxon contributes, so a per-taxon model matrix is not a row
## subset of the full-sample model matrix.
d4 = .make_data(seed = 4)
for (i in 1:20) d4$cnt[i, sample(seq_len(60), 12)] = 0
out$ancombc2_poly_term = .try(quote(
    ancombc2(data = d4$cnt, meta_data = d4$md,
             fix_formula = "poly(cont_cov, 2) + cat_cov", rand_formula = NULL,
             p_adj_method = "holm", pseudo_sens = FALSE,
             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
             group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
             alpha = 0.05, n_cl = 1, verbose = FALSE,
             global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)))

## Case 5: an unused level of the group factor reaches .get_struc_zero().
d5 = .make_data(seed = 5)
d5$md$cat_cov = factor(as.character(d5$md$cat_cov), levels = c("a", "b", "c", "z"))
out$struc_zero_unused_level = .try(quote(
    ANCOMBC:::.get_struc_zero(data = d5$cnt, meta_data = d5$md,
                              group = "cat_cov", neg_lb = FALSE)))
out$struc_zero_unused_level_neglb = .try(quote(
    ANCOMBC:::.get_struc_zero(data = d5$cnt, meta_data = d5$md,
                              group = "cat_cov", neg_lb = TRUE)))

## Case 6: a missing group label for some samples.
d6 = .make_data(seed = 6)
d6$md$cat_cov[c(1, 2, 31)] = NA
out$struc_zero_na_group = .try(quote(
    ANCOMBC:::.get_struc_zero(data = d6$cnt, meta_data = d6$md,
                              group = "cat_cov", neg_lb = TRUE)))

## Case 7: a two-level group, exercising the global test with a single
## contrasted coefficient.
d7 = .make_data(seed = 7)
d7$md$cat_cov = factor(ifelse(d7$md$cat_cov == "a", "a", "b"))
out$ancombc_two_level_group = .try(quote(
    ancombc(data = d7$cnt, meta_data = d7$md, formula = "cont_cov + cat_cov",
            pseudo = 1, p_adj_method = "holm", prv_cut = 0, lib_cut = 0,
            group = "cat_cov", struc_zero = TRUE, neg_lb = TRUE,
            tol = 1e-5, max_iter = 20, conserve = TRUE, alpha = 0.05,
            global = TRUE, n_cl = 1, verbose = FALSE)))

## Case 8: missing values in a covariate, so the design has incomplete rows.
d8 = .make_data(seed = 8)
d8$md$cont_cov[c(4, 17, 44)] = NA
for (i in 1:15) d8$cnt[i, sample(seq_len(60), 10)] = 0
out$ancombc2_na_covariate = .try(quote(
    ancombc2(data = d8$cnt, meta_data = d8$md,
             fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
             p_adj_method = "holm", pseudo_sens = FALSE,
             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
             group = "cat_cov", struc_zero = TRUE, neg_lb = TRUE,
             alpha = 0.05, n_cl = 1, verbose = FALSE,
             global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)))

for (nm in names(out)) {
    cat(sprintf("%-32s %s %s\n", nm, out[[nm]]$status,
                if (out[[nm]]$status == "error") out[[nm]]$msg else ""))
}
saveRDS(out, outfile)
cat("Wrote ", outfile, "\n", sep = "")

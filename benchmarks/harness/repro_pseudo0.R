# ---------------------------------------------------------------------------
# repro_pseudo0.R
#
# Usage: Rscript repro_pseudo0.R <libpath>
#
# Minimal reproduction of the ancombc(pseudo = 0) regression.
#
# ancombc() forms its response as y = log(feature_table + pseudo) and, unlike
# .ancombc2_core(), does not replace the resulting -Inf by NA. With pseudo = 0
# and a zero count, y therefore contains -Inf.
#
# Under 2.13.2 the per-taxon stats::lm() call is wrapped in try(), so the taxon
# with the non-finite response is returned as NA and the run completes.
# Under 2.15.1 the non-finite entry is treated as missing and the taxon is
# fitted on its observed samples.
#
# The script prints the error message if the call fails, and otherwise the log
# fold change of the affected taxon.
# ---------------------------------------------------------------------------

libpath = normalizePath(commandArgs(trailingOnly = TRUE)[1], mustWork = TRUE)
.libPaths(c(libpath, .libPaths()))
suppressPackageStartupMessages(library(ANCOMBC, lib.loc = libpath))
cat(sprintf("ANCOMBC %s\n", as.character(packageVersion("ANCOMBC"))))

set.seed(1)
cnt = matrix(rnbinom(20 * 30, mu = 200, size = 2), nrow = 20, ncol = 30)
rownames(cnt) = paste0("T", seq_len(20))
colnames(cnt) = paste0("S", seq_len(30))
cnt[1, 1:5] = 0                      # one taxon with zero counts
md = data.frame(grp = factor(rep(c("a", "b"), each = 15)),
                cont_cov = rnorm(30),
                row.names = colnames(cnt))

res = try(suppressWarnings(suppressMessages(
    ancombc(data = cnt, meta_data = md, formula = "grp + cont_cov", pseudo = 0,
            p_adj_method = "holm", prv_cut = 0, lib_cut = 0, group = "grp",
            struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, max_iter = 20,
            conserve = TRUE, alpha = 0.05, global = FALSE, n_cl = 1,
            verbose = FALSE))), silent = TRUE)

if (inherits(res, "try-error")) {
    cat("ERROR: ", conditionMessage(attr(res, "condition")), "\n", sep = "")
} else {
    cat("OK. lfc of the affected taxon T1:\n")
    print(res$res$lfc[res$res$lfc$taxon == "T1", ])
}

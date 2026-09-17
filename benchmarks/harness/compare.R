# ---------------------------------------------------------------------------
# compare.R
#
# Usage: Rscript compare.R <outdir_a> <outdir_b> [tolerance]
#
# For every workload with an RDS file in both directories, reports
# all.equal(a, b, tolerance = <tolerance>) (default 1e-8) and joins the timing
# and memory records from the two runs. The table is printed and written to
# <outdir_b>/compare_vs_<basename(outdir_a)>.csv. Exit status is 1 if any
# workload is not equal.
#
# Stripping: environment-bearing components are removed from both objects
# before comparison, because an environment carries the calling frame and
# compares unequal for reasons unrelated to the numerical result. The stripper
# drops any component that is a function, an environment, or carries an
# attribute environment (formula, terms, model fit objects), and every stripped
# path is printed. ANCOMBC's return values are lists of data frames, matrices
# and vectors, so in practice nothing is expected to be stripped; the pass is
# kept so that a future return value containing a fit object does not silently
# produce a spurious inequality.
# ---------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript compare.R <outdir_a> <outdir_b> [tolerance]")
}
dir_a = normalizePath(args[1], mustWork = TRUE)
dir_b = normalizePath(args[2], mustWork = TRUE)
tol = if (length(args) >= 3) as.numeric(args[3]) else 1e-8

stripped_paths = character(0)

.is_env_bearing = function(x) {
    is.function(x) || is.environment(x) ||
        !is.null(attr(x, ".Environment")) ||
        inherits(x, c("formula", "terms", "lm", "glm", "merMod", "lmerMod",
                      "call", "cluster"))
}

# Stripped components are replaced by a sentinel string rather than removed, so
# that list lengths and positions stay aligned between the two objects and a
# NULL component (an absent optional result) is preserved as NULL.
.strip = function(x, path = "") {
    if (.is_env_bearing(x)) {
        stripped_paths <<- c(stripped_paths,
                             sprintf("%s [%s]", path, class(x)[1]))
        return("<stripped: environment-bearing>")
    }
    if (is.list(x) && !is.data.frame(x)) {
        nms = names(x)
        for (i in seq_along(x)) {
            lbl = if (is.null(nms) || !nzchar(nms[i])) sprintf("[[%d]]", i) else
                paste0("$", nms[i])
            x[i] = list(.strip(x[[i]], paste0(path, lbl)))
        }
    }
    x
}

.read_timings = function(d) {
    f = file.path(d, "timings.csv")
    if (!file.exists(f)) return(NULL)
    read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

t_a = .read_timings(dir_a)
t_b = .read_timings(dir_b)

files_a = sub("\\.rds$", "", basename(Sys.glob(file.path(dir_a, "*.rds"))))
files_b = sub("\\.rds$", "", basename(Sys.glob(file.path(dir_b, "*.rds"))))
common = intersect(files_a, files_b)
only_a = setdiff(files_a, files_b)
only_b = setdiff(files_b, files_a)

if (length(common) == 0) stop("No workloads present in both directories.")
if (length(only_a) > 0)
    cat("Only in ", basename(dir_a), ": ", paste(only_a, collapse = ", "),
        "\n", sep = "")
if (length(only_b) > 0)
    cat("Only in ", basename(dir_b), ": ", paste(only_b, collapse = ", "),
        "\n", sep = "")

.lookup = function(tab, nm, col) {
    if (is.null(tab)) return(NA_real_)
    i = match(nm, tab$name)
    if (is.na(i)) return(NA_real_)
    as.numeric(tab[[col]][i])
}

rows = lapply(common, function(nm) {
    a = readRDS(file.path(dir_a, paste0(nm, ".rds")))
    b = readRDS(file.path(dir_b, paste0(nm, ".rds")))
    n0 = length(stripped_paths)
    a = .strip(a, nm)
    b = .strip(b, nm)
    if (length(stripped_paths) > n0) {
        cat("Stripped from ", nm, ": ",
            paste(unique(stripped_paths[(n0 + 1):length(stripped_paths)]),
                  collapse = "; "), "\n", sep = "")
    }
    cmp = all.equal(a, b, tolerance = tol)
    eq = isTRUE(cmp)
    msg = if (eq) "" else paste(as.character(cmp), collapse = " | ")
    ea = .lookup(t_a, nm, "elapsed_sec")
    eb = .lookup(t_b, nm, "elapsed_sec")
    data.frame(name = nm,
               equal = eq,
               first_diff = if (nchar(msg) > 300)
                   paste0(substr(msg, 1, 300), "...") else msg,
               elapsed_a = ea,
               elapsed_b = eb,
               speedup = round(ea / eb, 3),
               mem_a = .lookup(t_a, nm, "peak_mem_mb"),
               mem_b = .lookup(t_b, nm, "peak_mem_mb"),
               stringsAsFactors = FALSE)
})
out = do.call(rbind, rows)
out = out[order(out$name), ]

cat(sprintf("\nTolerance: %g\n\n", tol))
print(out, row.names = FALSE, right = FALSE)

csv = file.path(dir_b, paste0("compare_vs_", basename(dir_a), ".csv"))
write.csv(out, csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", csv))

n_bad = sum(!out$equal)
cat(sprintf("%d of %d workloads equal at tolerance %g\n",
            sum(out$equal), nrow(out), tol))
if (n_bad > 0) quit(status = 1)

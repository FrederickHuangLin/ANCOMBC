# ---------------------------------------------------------------------------
# run.R
#
# Usage: Rscript run.R <libpath> <outdir> [scale] [reps] [name_regex]
#
#   libpath : library directory holding the ANCOMBC installation to exercise.
#             ANCOMBC is attached from this path only.
#   outdir  : directory for <name>.rds and timings.csv. Created if absent.
#   scale   : "all" (default), "vignette" or "large".
#   reps    : number of timed repetitions per workload. Default 3 at vignette
#             scale and 1 at large scale. Pass "" or NA for the defaults.
#   name_regex : optional regular expression restricting the workload names.
#             When given, timings.csv is merged into any existing timings.csv in
#             outdir so that a partial re-run does not discard earlier rows.
#
# Wall-clock time is the median elapsed time of system.time() over reps.
# Peak memory is gc(reset = TRUE) before the call and gc() after, reporting the
# increase in "max used" Mb summed over Ncells and Vcells.
#
# NOTE: workloads that pass n_cl > 1 spawn worker processes (parallel PSOCK).
# gc() observes the master process only, so the reported peak memory excludes
# the workers. For those workloads the figure is a lower bound on the memory
# footprint of the run.
# ---------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript run.R <libpath> <outdir> [scale] [reps]")
}
libpath = normalizePath(args[1], mustWork = TRUE)
outdir = args[2]
scale_req = if (length(args) >= 3) args[3] else "all"
reps_req = if (length(args) >= 4 && nzchar(args[4])) as.integer(args[4]) else
    NA_integer_
name_regex = if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outdir = normalizePath(outdir, mustWork = TRUE)

# --- attach ANCOMBC from libpath only ---------------------------------------
# libpath is placed first so that the workloads' dependencies still resolve,
# and the loaded namespace is then verified to come from libpath.
.libPaths(c(libpath, .libPaths()))
suppressPackageStartupMessages(library(ANCOMBC, lib.loc = libpath))

pkg_version = as.character(utils::packageVersion("ANCOMBC"))
pkg_path = dirname(getNamespaceInfo("ANCOMBC", "path"))
if (!identical(normalizePath(pkg_path), libpath)) {
    stop(sprintf("ANCOMBC loaded from %s, expected %s", pkg_path, libpath))
}
cat(sprintf("ANCOMBC %s loaded from %s\n", pkg_version,
            getNamespaceInfo("ANCOMBC", "path")))
cat(sprintf("%s | host %s\n", R.version.string, Sys.info()[["nodename"]]))

script_dir = dirname(normalizePath(sub("^--file=", "", grep("^--file=",
                                   commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "workloads.R"))

sel = names(workloads)
if (!identical(scale_req, "all")) {
    sel = sel[vapply(workloads[sel], function(w) w$scale, character(1)) ==
                  scale_req]
}
if (!is.null(name_regex)) sel = grep(name_regex, sel, value = TRUE)
if (length(sel) == 0) stop("No workloads for scale '", scale_req, "'")

# --- fixed RNG state applied before every workload --------------------------
.reset_rng = function() {
    RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
            sample.kind = "Rejection")
    set.seed(123)
}

# Peak memory increase, in Mb, of the master process across a call, summed over
# Ncells and Vcells. The Mb figure is the unnamed column immediately after the
# "max used" cell count; its position varies with the platform (an R build with
# a memory limit inserts an extra column), so it is located by name.
.mem_col = function(g) {
    j = which(colnames(g) == "max used")
    if (length(j) != 1L) stop("cannot locate 'max used' column in gc() output")
    j + 1L
}
.mem_delta = function(g_before, g_after) {
    sum(g_after[, .mem_col(g_after)]) - sum(g_before[, .mem_col(g_before)])
}

rows = vector("list", length(sel))

for (k in seq_along(sel)) {
    nm = sel[k]
    w = workloads[[nm]]
    reps = if (!is.na(reps_req)) reps_req else if (w$scale == "large") 1L else 3L

    cat(sprintf("\n=== [%d/%d] %s (%s, %s, reps = %d)\n",
                k, length(sel), nm, w$scale, w$fn, reps))

    .reset_rng()
    d = w$prep()

    elapsed = numeric(reps)
    mem = numeric(reps)
    res = NULL
    for (r in seq_len(reps)) {
        .reset_rng()
        g0 = gc(reset = TRUE, full = TRUE)
        tm = system.time(res_r <- w$call(d))
        g1 = gc(full = TRUE)
        elapsed[r] = unname(tm[["elapsed"]])
        mem[r] = .mem_delta(g0, g1)
        if (r == 1L) res = res_r
        rm(res_r)
    }

    saveRDS(res, file.path(outdir, paste0(nm, ".rds")))
    rows[[k]] = data.frame(name = nm, scale = w$scale, `function` = w$fn,
                           elapsed_sec = round(median(elapsed), 3),
                           peak_mem_mb = round(median(mem), 1),
                           version = pkg_version,
                           check.names = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("--- %s: elapsed = %.2f s, peak mem = %.1f Mb\n",
                nm, median(elapsed), median(mem)))
    rm(d, res)
    gc(full = TRUE)
}

timings = do.call(rbind, rows)
csv = file.path(outdir, "timings.csv")
if (!is.null(name_regex) && file.exists(csv)) {
    old = read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    old = old[!(old$name %in% timings$name), , drop = FALSE]
    timings = rbind(old, timings)
    timings = timings[order(match(timings$name, names(workloads))), ]
}
write.csv(timings, csv, row.names = FALSE)
cat("\n")
print(timings, row.names = FALSE)
cat(sprintf("\nWrote %s\n", csv))

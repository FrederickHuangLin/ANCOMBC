# =============================================================================
# ANCOM-BC3 compiled computation engine
# =============================================================================
#
# ancombc3(engine = "compiled") evaluates the Poisson-log-normal quadrature of
# the second stage (.ac3_pln_obs(), ancombc3_pln.R) with the compiled routine
# ac3_pln_obs() of src/ancombc3_core.c instead of the R code. The routine
# performs the operations of the R code in the same order, with a table-based
# exponential whose relative error is of the order of one unit in the last
# place; the estimates and test statistics therefore agree with those of
# engine = "r" to within rounding error of that size propagated through the
# fit. engine = "r" uses the R code only.
#
# The selected engine is passed to the stage functions as their argument
# engine and stored in the per-taxon data of the taxon-wise fits
# (.ac3_da_taxon_data()), so that the inputs sent to worker processes carry
# it. The routine is part of the shared library of the package, which a
# worker process loads with the namespace of the package.

# Validated value of the argument engine.
.ac3_engine_check <- function(engine) {
  if (!is.character(engine) || length(engine) != 1L || !engine %in% c("r", "compiled")) {
    stop("engine must be \"r\" or \"compiled\".")
  }
  engine
}

# Compiled counterpart of .ac3_pln_obs(), with the same arguments and value.
# libm_exp = TRUE uses exp() of the C library for every exponential, which
# gives the values of the R code (used to check the routine).
.ac3_engine_pln_obs <- function(center, y, sigma, moments = TRUE, libm_exp = FALSE) {
  gl <- .ac3_pln_gl_nodes()
  o <- .Call(ac3_pln_obs_c, as.double(center), as.double(y), as.double(sigma),
             isTRUE(moments), gl$pos, gl$logw, .ac3_pln_level_nats, isTRUE(libm_exp))
  list(logf = o[[1L]], E1 = o[[2L]], E2 = o[[3L]], Ez = o[[4L]], Vz = o[[5L]])
}

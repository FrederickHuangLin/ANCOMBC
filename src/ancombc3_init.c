/* Interface between R and the core routines of ancombc3_core.c, and
 * registration of the routines for .Call(). */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#include "ancombc3_core.h"

SEXP ac3_mul_add_c(SEXP a_, SEXP b_, SEXP c_) {
  return ScalarReal(ac3_mul_add(asReal(a_), asReal(b_), asReal(c_)));
}

/* Poisson-log-normal cell quantities (.ac3_engine_pln_obs() in
 * R/ancombc3_engine.R). Returns the list (logf, E1, E2, Ez, Vz); the last
 * four elements are NULL when moments is FALSE. */
SEXP ac3_pln_obs_c(SEXP center_, SEXP y_, SEXP sigma_, SEXP moments_, SEXP pos_, SEXP logw_,
                   SEXP depth_, SEXP libm_exp_) {
  if (TYPEOF(center_) != REALSXP || TYPEOF(y_) != REALSXP || TYPEOF(pos_) != REALSXP ||
      TYPEOF(logw_) != REALSXP) {
    error("ac3_pln_obs_c: center, y, and the nodes must be double vectors.");
  }
  const R_xlen_t n = XLENGTH(center_);
  const int K = (int) XLENGTH(pos_);
  if (XLENGTH(y_) != n) error("ac3_pln_obs_c: center and y must have the same length.");
  if (XLENGTH(logw_) != K || K < 1) error("ac3_pln_obs_c: invalid quadrature nodes.");
  const double *center = REAL(center_), *y = REAL(y_), *pos = REAL(pos_), *logw = REAL(logw_);
  const double s = asReal(sigma_), depth = asReal(depth_);
  const int moments = asLogical(moments_), libm_exp = asLogical(libm_exp_);
  if (!R_FINITE(s) || s <= 0) error("ac3_pln_obs_c: sigma must be positive and finite.");
  if (moments == NA_LOGICAL || libm_exp == NA_LOGICAL) error("ac3_pln_obs_c: invalid logical argument.");
  for (R_xlen_t i = 0; i < n; i++) {
    if (!R_FINITE(center[i]) || !R_FINITE(y[i])) {
      error("ac3_pln_obs_c: center and y must be finite.");
    }
  }

  SEXP out = PROTECT(allocVector(VECSXP, 5));
  SEXP logf_ = allocVector(REALSXP, n);
  SET_VECTOR_ELT(out, 0, logf_);
  double *logf = REAL(logf_), *E1 = NULL, *E2 = NULL, *Ez = NULL, *Vz = NULL;
  if (moments) {
    SEXP e1_ = allocVector(REALSXP, n);
    SET_VECTOR_ELT(out, 1, e1_);
    SEXP e2_ = allocVector(REALSXP, n);
    SET_VECTOR_ELT(out, 2, e2_);
    SEXP ez_ = allocVector(REALSXP, n);
    SET_VECTOR_ELT(out, 3, ez_);
    SEXP vz_ = allocVector(REALSXP, n);
    SET_VECTOR_ELT(out, 4, vz_);
    E1 = REAL(e1_);
    E2 = REAL(e2_);
    Ez = REAL(ez_);
    Vz = REAL(vz_);
  }
  double *log_pois = (double *) R_alloc(n, sizeof(double));
  for (R_xlen_t i = 0; i < n; i++) log_pois[i] = dpois(y[i], y[i], 1);
  double *work = (double *) R_alloc(AC3_PLN_OBS_WORK(n, K), sizeof(double));
  ac3_pln_obs((ptrdiff_t) n, K, center, y, log_pois, s, moments, pos, logw, depth, libm_exp,
              logf, E1, E2, Ez, Vz, work);
  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef ac3_call_methods[] = {
  {"ac3_mul_add_c", (DL_FUNC) &ac3_mul_add_c, 3},
  {"ac3_pln_obs_c", (DL_FUNC) &ac3_pln_obs_c, 8},
  {NULL, NULL, 0}
};

void R_init_ANCOMBC(DllInfo *dll) {
  ac3_exp_table_init();
  R_registerRoutines(dll, NULL, ac3_call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

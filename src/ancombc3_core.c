/* Core routines of the ANCOM-BC3 compiled engine.
 *
 * This file uses plain C arguments and no R headers, so that the same source
 * is compiled by the R package and by the Python package. The interface to R
 * is in ancombc3_init.c.
 *
 * Floating-point expressions are evaluated as written, without contraction
 * into fused multiply-add instructions, so that the results do not depend on
 * the instruction set of the processor.
 *
 * ac3_pln_obs() evaluates the Poisson-log-normal cell quantities of
 * .ac3_pln_obs() in R/ancombc3_pln.R: log f(y; mu, sigma),
 * E(Lambda | y), E(Lambda^2 | y), E(z | y), and Var(z | y) for vectors center
 * and y of a common length and a scalar sigma, by the composite
 * Gauss-Legendre rule over the interval {v : L(v) >= L(v_hat) - depth} of the
 * log integrand in v = u - log A,
 *   L(v) = y {v - expm1(v)} - B expm1(v) - (v - d)^2 / (2 sigma^2),
 * with B = e^c for y = 0 and B = 0 for y > 0, A = y + B, d = c - log A,
 * Lambda = A e^v, and z = (v - d) / sigma.
 * The operations are those of the R code, in the same order:
 *   1. the mode v_hat of L by Newton-Raphson iterations applied to all cells
 *      until the largest sigma |L'(v)| is below 1e-11 (at most 60
 *      iterations), from the starting values of .ac3_pln_mode();
 *   2. the endpoints of the interval on each side of the mode by
 *      Newton-Raphson iterations applied to all cells until the largest
 *      absolute residual is below 1e-9 (at most 60 iterations), with each
 *      iterate kept strictly on its side of the mode;
 *   3. the log integrand at the nodes, its maximum over the nodes, and the
 *      sums over the nodes of the scaled integrand and of its products with
 *      Lambda, Lambda^2, z, and (z - E(z | y))^2, accumulated in long double
 *      as in rowSums().
 * The exponential is computed by ac3_exp(): argument reduction
 * x = (k / 128) log 2 + r with |r| <= log(2) / 256, a table of 2^(j / 128),
 * j = 0, ..., 127, and the degree-five Taylor polynomial of e^r - 1, whose
 * truncation error is below 1e-18 relative; the result is within about one
 * unit in the last place of e^x. Arguments outside (-700, 700) and NaN are
 * passed to exp() of the C library. expm1(v) is expm1() of the C library for
 * |v| < 0.25 and ac3_exp(v) - 1 otherwise. With libm_exp = TRUE every
 * exponential is exp() or expm1() of the C library, and the results equal
 * those of the R code.
 *
 * ac3_exp_table_init() fills the table of the exponential and is called once
 * before the first call of ac3_pln_obs(). */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract (off)
#endif

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "ancombc3_core.h"

#define AC3_LOG_SQRT_2PI 0.918938533204672741780329736406
#define AC3_TABLE_SIZE 128
#define AC3_MAX_ITER 60

static double ac3_exp_table[AC3_TABLE_SIZE];

void ac3_exp_table_init(void) {
  for (int j = 0; j < AC3_TABLE_SIZE; j++) {
    ac3_exp_table[j] = pow(2.0, (double) j / AC3_TABLE_SIZE);
  }
}

/* e^x by table-based argument reduction. k = round(128 x / log 2) is formed
 * by adding and subtracting 1.5 * 2^52; log(2) / 128 is split into a high part
 * with 36 significant bits, whose product with k is exact for |k| < 2^17, and
 * a low part. With k = 128 e + j, 0 <= j < 128,
 * e^x = 2^e 2^(j / 128) e^r, and 2^e is formed from its exponent bits. */
static inline double ac3_exp_table_based(double x) {
  static const double inv_step = 0x1.71547652b82fep+0 * AC3_TABLE_SIZE;
  static const double step_hi = 0x1.62e42fefa0000p-8;
  static const double step_lo = 0x1.cf79abc9e3b3ap-47;
  if (!(x > -700.0 && x < 700.0)) return exp(x);
  double kd = (x * inv_step + 0x1.8p52) - 0x1.8p52;
  int64_t k = (int64_t) kd;
  double r = x - kd * step_hi - kd * step_lo;
  double r2 = r * r;
  double p = r + r2 * (0.5 + r * (1.0 / 6.0)) + r2 * r2 * (1.0 / 24.0 + r * (1.0 / 120.0));
  int64_t j = k & (AC3_TABLE_SIZE - 1);
  int64_t e = (k - j) / AC3_TABLE_SIZE;
  double t = ac3_exp_table[j];
  double v = t + t * p;
  uint64_t bits = (uint64_t) (e + 1023) << 52;
  double scale;
  memcpy(&scale, &bits, sizeof(double));
  return v * scale;
}

static inline double ac3_exp(double x, int libm_exp) {
  return libm_exp ? exp(x) : ac3_exp_table_based(x);
}

static inline double ac3_expm1(double v, int libm_exp) {
  return (libm_exp || fabs(v) < 0.25) ? expm1(v) : ac3_exp_table_based(v) - 1.0;
}

/* L'(v) = -y expm1(v) - B e^v - (v - d) / sigma^2 (.ac3_pln_score_v()). */
static inline double ac3_score_v(double v, double em1, double y, double B, double d, double s2) {
  return -y * em1 - B * (1 + em1) - (v - d) / s2;
}

/* L(v) without its constant (.ac3_pln_L()). */
static inline double ac3_L(double v, double em1, double y, double B, double d, double s) {
  double zq = (v - d) / s;
  return y * (v - em1) - B * em1 - 0.5 * (zq * zq);
}

/* Returns a * b + c with the product rounded before the addition. Used to
 * verify that the compiler does not contract the expression. */
double ac3_mul_add(double a, double b, double c) {
  return a * b + c;
}

/* center, y, and log_pois have length n; pos and logw have length K. log_pois
 * holds the log Poisson mass at y with mean y (dpois(y, y, log = TRUE)). The
 * results logf and, with moments != 0, E1, E2, Ez, and Vz have length n; with
 * moments == 0 the last four are not accessed. work has length
 * AC3_PLN_OBS_WORK(n, K). The arguments are validated by the caller. */
void ac3_pln_obs(ptrdiff_t n, int K, const double *center, const double *y, const double *log_pois,
                 double s, int moments, const double *pos, const double *logw, double depth,
                 int libm_exp, double *logf, double *E1, double *E2, double *Ez, double *Vz,
                 double *work) {
  double *lg = work;
  double *lamk = lg + K;
  double *zk = lamk + K;
  double *ek = zk + K;
  double *Bv = ek + K;
  double *dv = Bv + n;
  double *cst = dv + n;
  double *vhat = cst + n;
  double *curv = vhat + n;
  double *Lhat = curv + n;
  double *v_lo = Lhat + n;
  double *v_hi = v_lo + n;
  const double s2 = s * s, log_s = log(s);

  /* Reference intensity (.ac3_pln_reference()) and starting value of the mode
   * (.ac3_pln_mode()). */
  for (ptrdiff_t i = 0; i < n; i++) {
    if (y[i] > 0) {
      Bv[i] = 0.0;
      dv[i] = center[i] - log(y[i] > 1e-300 ? y[i] : 1e-300);
      double shifted = dv[i] + s2 * y[i];
      vhat[i] = dv[i] < 0 ? (shifted < 0 ? shifted : 0.0) : dv[i];
    } else {
      Bv[i] = exp(center[i]);
      dv[i] = 0.0;
      vhat[i] = -log1p(s2 * Bv[i]);
    }
    cst[i] = log_pois[i] - Bv[i];
  }

  /* 1. Mode of L(v). curv and Lhat hold expm1(v) and the score during the
   * iterations. */
  for (int it = 0; it < AC3_MAX_ITER; it++) {
    double largest = 0;
    for (ptrdiff_t i = 0; i < n; i++) {
      double em1 = ac3_expm1(vhat[i], libm_exp);
      double g = ac3_score_v(vhat[i], em1, y[i], Bv[i], dv[i], s2);
      curv[i] = em1;
      Lhat[i] = g;
      if (fabs(s * g) > largest) largest = fabs(s * g);
    }
    if (largest < 1e-11) break;
    for (ptrdiff_t i = 0; i < n; i++) {
      vhat[i] = vhat[i] + Lhat[i] / ((y[i] + Bv[i]) * (1 + curv[i]) + 1 / s2);
    }
  }
  for (ptrdiff_t i = 0; i < n; i++) {
    double em1 = ac3_expm1(vhat[i], libm_exp);
    curv[i] = s2 * (y[i] + Bv[i]) * (1 + em1) + 1;
    Lhat[i] = ac3_L(vhat[i], em1, y[i], Bv[i], dv[i], s);
  }

  /* 2. Endpoints of the integration interval (.ac3_pln_level_boundary()). */
  for (int side_index = 0; side_index < 2; side_index++) {
    const double side = side_index == 0 ? -1.0 : 1.0;
    double *v = side_index == 0 ? v_lo : v_hi;
    for (ptrdiff_t i = 0; i < n; i++) v[i] = vhat[i] + side * s * sqrt(2 * depth / curv[i]);
    for (int it = 0; it < AC3_MAX_ITER; it++) {
      double largest = 0;
      for (ptrdiff_t i = 0; i < n; i++) {
        double em1 = ac3_expm1(v[i], libm_exp);
        double h = ac3_L(v[i], em1, y[i], Bv[i], dv[i], s) - (Lhat[i] - depth);
        if (fabs(h) > largest) largest = fabs(h);
      }
      if (largest < 1e-9) break;
      for (ptrdiff_t i = 0; i < n; i++) {
        double em1 = ac3_expm1(v[i], libm_exp);
        double h = ac3_L(v[i], em1, y[i], Bv[i], dv[i], s) - (Lhat[i] - depth);
        double h1 = ac3_score_v(v[i], em1, y[i], Bv[i], dv[i], s2);
        double v_new = v[i] - h / h1;
        double bound = side < 0 ? vhat[i] - s * 1e-6 / sqrt(curv[i]) : vhat[i] + s * 1e-6 / sqrt(curv[i]);
        if (side < 0) {
          v[i] = v_new < bound ? v_new : bound;
        } else {
          v[i] = v_new > bound ? v_new : bound;
        }
      }
    }
  }

  /* 3. Quadrature in v over [v_lo, v_lo + width]. */
  for (ptrdiff_t i = 0; i < n; i++) {
    const double yi = y[i], Bi = Bv[i], di = dv[i], Ai = y[i] + Bv[i];
    const double width = v_hi[i] - v_lo[i];
    const double log_width = log(width);
    double m = -INFINITY;
    for (int k = 0; k < K; k++) {
      double v = v_lo[i] + width * pos[k];
      double em1 = ac3_expm1(v, libm_exp);
      double zq = (v - di) / s;
      double val = yi * (v - em1) - Bi * em1 + cst[i] + (-(AC3_LOG_SQRT_2PI + 0.5 * zq * zq)) - log_s +
        log_width + logw[k];
      lg[k] = val;
      lamk[k] = Ai * (1 + em1);
      zk[k] = zq;
      if (val > m) m = val;
    }
    long double a0 = 0, a1 = 0, a2 = 0, az = 0, avz = 0;
    for (int k = 0; k < K; k++) {
      double e = ac3_exp(lg[k] - m, libm_exp);
      ek[k] = e;
      a0 += e;
      if (moments) {
        a1 += e * lamk[k];
        a2 += e * (lamk[k] * lamk[k]);
        az += e * zk[k];
      }
    }
    double s0 = (double) a0;
    logf[i] = m + log(s0);
    if (moments) {
      E1[i] = (double) a1 / s0;
      E2[i] = (double) a2 / s0;
      Ez[i] = (double) az / s0;
      for (int k = 0; k < K; k++) {
        double dz = zk[k] - Ez[i];
        avz += ek[k] * (dz * dz);
      }
      Vz[i] = (double) avz / s0;
    }
  }
}

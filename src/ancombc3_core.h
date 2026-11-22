#ifndef ANCOMBC3_CORE_H
#define ANCOMBC3_CORE_H

#include <stddef.h>

/* Length of the work array of ac3_pln_obs() for n cells and K nodes. */
#define AC3_PLN_OBS_WORK(n, K) (8 * (size_t) (n) + 4 * (size_t) (K))

void ac3_exp_table_init(void);

double ac3_mul_add(double a, double b, double c);

void ac3_pln_obs(ptrdiff_t n, int K, const double *center, const double *y, const double *log_pois,
                 double s, int moments, const double *pos, const double *logw, double depth,
                 int libm_exp, double *logf, double *E1, double *E2, double *Ez, double *Vz,
                 double *work);

#endif

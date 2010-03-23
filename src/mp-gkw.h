/*
 * mp-gkw.h
 *
 * FUNCTION:
 * Compute matrix elts of the GKW operator.
 *
 * HISTORY:
 * Linas Jan 2010
 */

#include <gmp.h>

#ifdef  __cplusplus
extern "C" {
#endif

// Return the matrix element for the matrix element G_mp of the GKW
// operator, expanded at the x=1 location.
void gkw(mpf_t elt, int m, int p, unsigned int prec);

// Return the continuous-valued version of the GKW operator.
// (the matrix elts occur at integer values)
// This implementation uses GMP multi-precision
void gkw_smooth(mpf_t elt, double m, double p, unsigned int prec);

#ifdef  __cplusplus
};
#endif


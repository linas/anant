/*
 * Generating functions for assorted number-theoretic functions
 * Implementation in bignums.
 *
 * October 2016
 */

#include <mp-complex.h>

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * Ordinary generating function for arithmetic series.
 */
void cpx_ordinary_genfunc(cpx_t sum, cpx_t z, int prec,
                          long (*func)(long));

/*
 * Exponential generating function for arithmetic series.
 * The second form expects func to return a value in the reference
 * mpf_t*
 */
void cpx_exponential_genfunc(cpx_t sum, cpx_t z, int prec,
                             long (*func)(long));
void cpx_exponential_genfunc_mpf(cpx_t sum, cpx_t z, int prec,
                                 void (*func)(mpf_t, long));

/*
 * Dyadic Farey rescaling.
 */
void cpx_exponential_twist(cpx_t sum, cpx_t z, int prec,
                           long (*func)(long));
#ifdef  __cplusplus
};
#endif

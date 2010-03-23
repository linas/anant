/*
 * mp-binomial.h
 *
 * High-precison factorials and binomials, using the 
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2005,2006 Linas Vepstas
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

#include <gmp.h>
#include "mp-complex.h"

#ifdef  __cplusplus
extern "C" {
#endif

/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */
void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n);

/** 
 * i_factorial -- the factorial
 */
// #define USE_LOCAL_FACTORIAL
#ifdef USE_LOCAL_FACTORIAL
void i_factorial (mpz_t fact, unsigned int n);
#else
#define i_factorial mpz_fac_ui
#endif /* USE_LOCAL_FACTORIAL */

/** 
 * fp_inv_factorial -- Return 1/n!
 * Returns a cached value to improve performance
 */
void fp_inv_factorial (mpf_t invfac, unsigned int n, unsigned int prec);

/* i_binomial
 * Binomial coefficient (n k)
 */
// #define USE_LOCAL_BINOMIAL
#ifdef USE_LOCAL_BINOMIAL
void i_binomial (mpz_t bin, unsigned int n, unsigned int k);
#else 
#define i_binomial mpz_bin_uiui
#endif /* USE_LOCAL_BINOMIAL */

/**
 * i_binomial_sequence -- returns binomial, assumes purely sequential access
 * 
 * This routine assumes that the binomial coefficients will be 
 * accessed in an utterly sequential mode, with k running from 
 * zero to n, and n running from zero to k. For sequential access,
 * this routine is very very fast. Otherwise, random access is used
 * which is considerably slower.
 */
void i_binomial_sequence (mpz_t bin, unsigned int n, unsigned int k);

/** 
 * stirling_first - Stirling Numbers of the First kind, 
 * normalized so that they are all positive.
 * Uses dynamically-sized cache.
 */
void i_stirling_first (mpz_t s, unsigned int n, unsigned int k);
/* A funny off-by-one sum of stirling and binomial */
void i_stirbin_sum (mpz_t s, unsigned int n, unsigned int m);

/** 
 * stirling_second - Stirling Numbers of the Second kind, 
 * Uses dynamically-sized cache.
 */
void i_stirling_second (mpz_t s, unsigned int n, unsigned int k);

/* binomial transform of power sum */
void fp_bin_xform_pow (mpf_t bxp, unsigned int n, unsigned int s);

/** 
 * fp_harmonic -- The harmonic number
 */
void fp_harmonic (mpf_t harm, unsigned int n, unsigned int prec);

/** 
 * fp_poch_rising
 * Rising pochhammer symbol (x)_n, for real values of x and integer n.
 * p = x(x+1)(x+2)...(x+n-1) = Gamma(x+n)/Gamma(x)
 *
 * Brute force, simple.
 */
void fp_poch_rising (mpf_t poch, mpf_t x, unsigned int n);
void fp_poch_rising_d (mpf_t poch, double x, unsigned int n);

/* cpx_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */
void cpx_poch_rising_d (cpx_t poch, double re_s, double im_s, unsigned int n);

/* cpx_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */
void cpx_poch_rising (cpx_t poch, const cpx_t s, unsigned int n);

/* fp_binomial
 * Binomial coefficient 
 */
void fp_binomial_d (mpf_t bin, double s, unsigned int k);

/**
 * cpx_binomial-- Complex binomial coefficient
 * Compute the binomial coefficient (s, k)
 */
void cpx_binomial_d (cpx_t bin, double re_s, double im_s, unsigned int k);
void cpx_binomial (cpx_t bin, const cpx_t s, unsigned int k);

/**
 * cpx_binomial_sum_cache-- Complex binomial coefficient
 * Compute the binomial coefficient for (s+k, k)
 * This routine caches computed values, and thus can be considerably
 * faster if called again with the same k and s values.
 */
void cpx_binomial_sum_cache (cpx_t bin, const cpx_t s, unsigned int k);

#ifdef  __cplusplus
};
#endif


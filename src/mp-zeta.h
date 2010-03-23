/*
 * mp-zeta.h
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Copyright (C) 2005 Linas Vepstas
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

/* Fixed-point bernoulli number */
void q_bernoulli (mpq_t bern, int n);

/* Compute and return the "exact" result for the zeta function for 
 * any value of even n. Computed top prec decimal places.
 * works by computing the Bernoulli number first.
 */
void fp_zeta_even (mpf_t zeta, unsigned int n, int prec);

/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Carries out math to prec decimal digits
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec);

/* Same, using Helmut Hasse convergent algo */
void fp_hasse_zeta (mpf_t zeta, unsigned int s, int prec);

/* Same, using P. Borwein convergent algo */
void fp_borwein_zeta (mpf_t zeta, unsigned int s, int prec);

/**
 * cpx_borwein_zeta -- use the Borwein algorithm for complex argument
 * 
 * Compute and return a value for the Riemann zeta function, for
 * a complex value of 's'. Uses the  P. Borwein algorithm for 
 * rapid computation.
 */
void cpx_borwein_zeta (cpx_t zeta, const cpx_t ess, int prec);

/**
 * cpx_borwein_zeta_cache -- Caching Riemann zeta for complex argument
 * 
 * Compute and return a value for the Riemann zeta function, for
 * a complex value of 's+n'. Uses the  P. Borwein algorithm for 
 * rapid computation.
 *
 * If the value of 's' is held constant, while the value of 'n'
 * is varied, then this routine caches the computed values,
 * returning the cached values on the second and later calls,
 * thus avoiding the overhead of repeated recalculation.
 */
void cpx_borwein_zeta_cache (cpx_t zeta, const cpx_t ess, unsigned int n, int prec);

/* Brute-force summation */
void fp_zeta_brute (mpf_t zeta, unsigned int s, int prec);

/* Stieltjes constants */
// void stieltjes_gamma (mpf_t gam, int n);

/* 
 * Compute a_sub_n
 * the w argument is for the power bit -- 
 */
void a_sub_n (mpf_t a_n, mpf_t w, unsigned int n, unsigned int prec);
void b_sub_n (mpf_t b_n, unsigned int n, unsigned int prec);

/* compute a_sub_s for complex-valued s
 */
void a_sub_s (mpf_t re_a, mpf_t im_a, double re_s, double im_s, unsigned int prec);
void b_sub_s_d (mpf_t re_a, mpf_t im_a, double re_s, double im_s, 
              unsigned int prec, int nterms, double eps);
void b_sub_s (mpf_t re_a, mpf_t im_a, mpf_t re_s, mpf_t im_s, 
              unsigned int prec, int nterms, double eps);

#ifdef  __cplusplus
};
#endif


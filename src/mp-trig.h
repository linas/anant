/*
 * mp-trig.h
 *
 * High-precison Elementary functions, using the 
 * Gnu Multiple-precision library.
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
 *
 */

#include <gmp.h>
#include "mp-complex.h"

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * i_pow - raise n to the m power
 */

void i_pow (mpz_t p, unsigned int n, unsigned int m);

/**
 * fp_inv_pow - raise n to the -m power, where m must be positive. 
 */
void fp_inv_pow (mpf_t p, unsigned int n, unsigned int m);

/**
 * fp_exp -  Floating point exponential
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants. 
 *
 * The complex exp is built up from the real trig functions.
 * The complex trig functions are built up from the complex exp.
 * In all cases, the most basic idents are used, so these
 * are not speedy!.
 */
void fp_exp (mpf_t ex, const mpf_t z, unsigned int prec);
void fp_sine (mpf_t sine, const mpf_t z, unsigned int prec);
void fp_cosine (mpf_t cosine, const mpf_t z, unsigned int prec);
void cpx_exp (cpx_t ex, const cpx_t z, unsigned int prec);
void cpx_sine (cpx_t sine, const cpx_t z, unsigned int prec);
void cpx_cosine (cpx_t cosine, const cpx_t z, unsigned int prec);
void cpx_tangent (cpx_t tang, const cpx_t z, unsigned int prec);

/**
 * fp_log -  Floating point logarithm
 * Implemented using a brute-force, simple algo, with 
 * minor attempts at optimization. 
 *
 * fp_log_m1 computes -log(1-z) using Taylor's expansion for small z.
 * Does not perform any other optimizations -- just simply sums the
 * Taylor series, and that's all.
 *
 * fp_log_ui takes integer arguments, and keeps previous 
 * results cached for improved performance.
 */
void fp_log_m1 (mpf_t lg, const mpf_t z, unsigned int prec);
void fp_log (mpf_t lg, const mpf_t z, unsigned int prec);
void fp_log_ui (mpf_t lg, unsigned int z, unsigned int prec);
void cpx_log_m1 (cpx_t lg, const cpx_t z, unsigned int prec);
void cpx_log (cpx_t lg, const cpx_t z, unsigned int prec);

/**
 * fp_arctan -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Very slow near y=x
 */
void fp_arctan (mpf_t atn, const mpf_t z, unsigned int prec);
void fp_arctan2 (mpf_t atn, const mpf_t y, const mpf_t x, unsigned int prec);

/**
 * cpx_sqrt
 * Simple implemenmtation of complex square-root
 */
void cpx_sqrt (cpx_t sqrt, const cpx_t zee, int prec);

/**
 * cpx_pow-- return q^s for complex q, s.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ...
 */
void cpx_pow (cpx_t powc, const cpx_t q, const cpx_t ess, int prec);

/**
 * cpx_mpf_pow-- return q^s for complex s, real, positive q.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ...
 */
void cpx_mpf_pow (cpx_t powc, const mpf_t q, const cpx_t ess, int prec);

/**
 * cpx_pow_ui-- return q^n for complex q, positive integer n.
 * Should be fairly speedy, as it uses a log(n) implementation
 */
void cpx_pow_ui (cpx_t powc, const cpx_t q, unsigned int n);

/**
 * cpx_ui_pow -- return k^s for complex s, integer k.
 *
 * Uses a brute-force algo: it requires a logarithm, an exp, sin 
 * and cos to be computed, each of which are kinda slow ... 
 */
void cpx_ui_pow (cpx_t powc, unsigned int k, const cpx_t ess, int prec);

/**
 * cpx_ui_pow_cache -- return k^s for complex s, integer k.
 *
 * If s is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds. If is is changed from call
 * to call, then the cache is cleared.
 */
void cpx_ui_pow_cache (cpx_t powc, unsigned int k, const cpx_t ess, int prec);

/**
 * fp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */

void fp_pow_rc (cpx_t diri, int k, const mpf_t q, const cpx_t ess, int prec);
void cpx_pow_rc (cpx_t diri, int k, const cpx_t q, const cpx_t ess, int prec);

#ifdef  __cplusplus
};
#endif

/* =============================== END OF FILE =========================== */


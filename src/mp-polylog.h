/** 
 * mp-polylog.h
 *
 * Implement Borwein-style polylogarithm.
 * Also implement the "periodic zeta" and 
 * the Hurwitz zeta function.
 *
 * Copyright (C) 2006,2007 Linas Vepstas
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
#include "mp-complex.h"

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * cpx_polylog_nint -- compute the polylogarithm at negetive integers
 *
 * Li_{-n}(z) 
 * At the negative integers, the polylog is a rational function,
 * meromorphic everywhere except for multiple poles at z=1.
 */
void cpx_polylog_nint (cpx_t plog, unsigned int negn, const cpx_t zee);

/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * The magnitude of z must be less than one in order for the 
 * summation to be caqrried out.
 *
 * Caches intermediate results, so that overall performance is
 * considerably better if z is varied while s is held fixed.
 */
void cpx_polylog_sum (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec);

/**
 * cpx_polylog -- polylogarithm
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * Works for general complex s, z; lightly tested, may be buggy.
 * Watch out for branchpoint at z=1.
 *
 * Returns a non-zero value if algo was unable to evaluate at
 * the given point.
 */
int cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec);

/**
 * cpx_polylog_euler -- compute the polylogarithm from Hurwitz Euler.
 *
 * Combine two Hurwitz Euler-Maclaurin evaluations to obtain the polylogarithm.
 */
void cpx_polylog_euler (cpx_t zeta, const cpx_t ess, const cpx_t zee, int prec);

/**
 * cpx_polylog_sheet -- give the branch difference for the polylog
 * M is the monodromy number of going around z=0
 * N is the monodromy number of going around z=1
 *
 * For M=0, the branch difference is that given between the principle
 * sheet, and the N'th winding around the z=1 branch point.
 * Thus, for example, the (0,1)'th sheet of Li_s(z) is given by
 *      (2pi i)^s (ln z/(2pi i)^{s-1} / Gamma (s)
 *
 * For N=0, the monodromy M has no effect.
 * For N!=0, the monodromy is given with respect to the N=1 sheet.
 */
void cpx_polylog_sheet(cpx_t delta, const cpx_t ess, const cpx_t zee, int z0_dromy, int z1_dromy, int prec);
void cpx_polylog_sheet_g0_action(cpx_t delta, const cpx_t ess, int direction, int prec);
void cpx_polylog_sheet_g1_action(cpx_t delta, const cpx_t ess, const cpx_t zee, int sheet, int direction, int prec);


/**
 * cpx_periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 */
void cpx_periodic_zeta (cpx_t z, const cpx_t ess, const mpf_t que, int prec);

/**
 * cpx_periodic_beta -- Periodic beta function 
 *
 * Similar to periodic zeta, but with different normalization
 *
 * beta = 2 Gamma(s+1) (2\pi)^{-s} F(s,q)
 *
 * Caches intermediate terms, and so performance is much better 
 * if s is held const, while q is varied.
 */
void cpx_periodic_beta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec);

/**
 * cpx_hurwitz_zeta -- Hurwitz zeta function
 * Returns zeta = sum_{n=0}^infty 1/(n+q)^s
 * Accepts complex s, real-valued q.
 *
 * Built up from the fast polylogarithm algo
 * Caches intermediate terms, and so performance is much better 
 * if s is held const, while q is varied.
 * The input value of q must be postive; 
 * the algo gets slow if q is very large.
 */
void cpx_hurwitz_zeta (cpx_t hzeta, const cpx_t ess, const mpf_t que, int prec);

/**
 * cpx_hurwitz_taylor -- Hurwitz zeta function taylor series
 *
 * Implement the Hurwitz zeta as a taylor expansion about q=0
 * (pulling out the leading 1/q^s term to handle uniquely)
 */
void cpx_hurwitz_taylor (cpx_t hzeta, const cpx_t ess, const cpx_t que, int prec);

/**
 * cpx_hurwitz_euler -- Hurwitz zeta function via Euler-Maclaurin algo
 *
 * This function computes the value of the Hurwitz zeta function
 * using an Euler-Maclaurin summation to obtain an estimate.
 *
 * The algorithm appears to work in principle (well, it gets 4 or 5
 * digits right), but we are doing a really really bad error estimate.
 * So it doesn't work at higher precision, at least not on the critical
 * strip.
 */
void cpx_hurwitz_euler_fp(cpx_t hzeta, cpx_t ess, mpf_t que, int prec);
void cpx_hurwitz_euler(cpx_t hzeta, cpx_t ess, cpx_t que, int prec);

#ifdef  __cplusplus
};
#endif


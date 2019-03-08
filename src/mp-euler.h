/**
 * mp-euler.h
 *
 * Euler resummation of slowly convergent (alternating) sequences.
 *
 * Copyright (C) 2019 Linas Vepstas
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

#ifndef __MP_EULER_H__
#define __MP_EULER_H__

#ifdef  __cplusplus
extern "C" {
#endif

#include "mp-complex.h"

/* =============================================== */
/**
 * cpx_euler_sum
 *
 * Implement Euler (re-)summation (aka Euler transformation) for
 * complex-valued arithmetic functions. That is, given f(n) defined
 * on positive integers n, return
 *
 * $ sum_{n=0}^\infty 2^{-(n+1)} \sum_{k=0}^n {n \choose k} f(k+1) $
 *
 * Note the sign convention used here: no minus signs; so that sign
 * alternation (if any) is in the f(n).
 *
 * See Wikipedia for more about Euler summation.
 *
 * See also Jonathan Sondow, "Analytic continuation of Riemann's zeta
 * function and values at negative integers via Euler's transformation
 * of series" (1994) Proceedings of the American Mathematical Society,
 * Vol. 120, pp 421-424.
 * http://www.ams.org/journals/proc/1994-120-02/S0002-9939-1994-1172954-7/S0002-9939-1994-1172954-7.pdf
 * for a particularly simple and direct statement, including a direct
 * application.
 *
 * The Euler transformation is an extremely simple kind of resummation.
 * In certain cases, it can work marvels for conditionally convergent
 * alternating series. Far more powerful resummation techniques exist;
 * see, for example, the references to P. Borwein convergent algo given
 * in the polylog code.
 *
 * The Euler transformation can also completely fail to accelerate
 * series convergence, if the series is not quite of the right form.
 * I am not aware of any results indicating when it works, and when it
 * fails, or why. Its a fairly serious question, having reprecussions
 * in number theory. (Viz, why the Hasse/Sondow resummation converges
 * rapidly and uniformly, while nearby sequences completely fail to do
 * so.)
 *
 * @func function providing values to be resummed.
 *       It will only be called for positive integers.
 *       'nprec' is the suggested decimal precisiton at which 'fun'
 *       should perform its calculations.
 * @ndigits Terminate summation when each term contributes less than
 *       10^ndigits to the sum. This is approximately that accuracy
 *       that will be attained. Roughly.
 * @maxterms Maximum number of terms to sum. If few terms are needed
 *       to reach the desired precision, that's fine; summation will
 *       stop. Otherwise, is summation has not converged, it will be
 *       terminated at `maxterms`.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 *
 * @returns number of terms actually summed.
 */
unsigned int cpx_euler_sum(cpx_t result,
              void (*func)(cpx_t f, unsigned long p, int nprec),
              unsigned int ndigits,
              unsigned int maxterms,
              int nprec);

/* =============================================== */
/**
 * cpx_newton_series
 *
 * Implement a Newton series interpolation for complex-valued arithmetic
 * functions. That is, given f(n) defined on positive integers n, return
 *
 * $ f(z) = sum_{n=0}^\infty (-1)^{-n} {z-1 \choose n}
 *       \sum_{k=0}^n (-1)^{-k} {n \choose k} f(k+1) $
 *
 * Note the sign convention above: there are two sources of alternating
 * signs.
 *
 * See Wikipedia for more about finite differences and Newton series.
 *
 * See also Philippe Flajolet and Robert Sedgewick, "Mellin Transforms
 * and Asymptotics: Finite Differences and Rice's Integrals" (1995)
 * Theoretical Computer Science, vol. 144 pages 101-124.
 * https://www.sciencedirect.com/science/article/pii/030439759400281M
 * for general thoughts.
 *
 */
unsigned int cpx_newton_series(cpx_t result,
              void (*func)(cpx_t f, unsigned long p, int nprec),
              cpx_t zee,
              unsigned int ndigits,
              unsigned int maxterms,
              int nprec);

#ifdef  __cplusplus
};
#endif

#endif /* __MP_EULER_H__ */

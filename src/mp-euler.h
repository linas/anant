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
 * cpx_euler.
 *
 * Implement Euler summation for complex-valued arithmetic function.
 * That is, given f(n) defined on positive integers n, return
 *
 * $ sum_{n=0}^\infty 2^{-(n+1)} \sum_{k=0}^n {n \choose k} f(k+1) $
 *
 * See Wikipedia for more about Euler summation.
 *
 * This is an extremely simple kind of resummation. It usually works
 * marvels for conditionally convergent alternating series.
 * Far more powerful resummation techniques exist; see, for example,
 * the references to P. Borwein convergent algo given in the polylog
 * code.
 *
 * @func function providing values to be resummed.
 *       It will only be called for positive integers.
 *       'nprec' is the suggested decimal precisiton at which 'fun'
 *       should perform its calculations.
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
              unsigned int maxterms,
              int nprec);

#ifdef  __cplusplus
};
#endif

#endif /* __MP_EULER_H__ */

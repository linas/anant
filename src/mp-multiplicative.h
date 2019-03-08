/**
 * mp-multiplicative.h
 *
 * Fill in values of completely multiplicative arithmetic function.
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

#ifndef __MP_MULTIPLICATIVE_H__
#define __MP_MULTIPLICATIVE_H__

#ifdef  __cplusplus
extern "C" {
#endif

#include "mp-complex.h"

/* =============================================== */
/**
 * cpx_multiplicative.
 *
 * Fill in values of completely multiplicative arithmetic function.
 * That is, given a function defined only on the prime numbers,
 * provide a values for on composite integers, by factoring them.
 *
 * A completely multiplicative function is a complex-valued function
 * on the positive integers that is a homomorpism under multiplication
 * i.e. preserves multiplication, i.e. f(mn) = f(m) f(n) for positive
 * integers m,n. Thus, it is sufficient to specify f(p) for prime p,
 * all other values are determined by the values on the primes.
 *
 * Given only f(p), this computes the values at all other integers.
 *
 * Two variants are provided: a cachine and a non-caching version.
 * The caching version will call f(p) only one per prime p. The
 * non-caching version might call f(p) repeatedly.
 *
 * @func function provding values on the prime numbers.
 *       It will only be called for a prime-number argument.
 *       'nprec' is the suggested decimal precisiton at which 'fun'
 *       should perform its calculations.
 * @n    number at which results should be computed.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 */
void cpx_multiplicative(cpx_t result,
              void (*func)(cpx_t f, unsigned long p, int nprec),
              unsigned long n,
              int nprec);

#ifdef  __cplusplus
};
#endif

#endif /* __MP_MULTIPLICATIVE_H__ */

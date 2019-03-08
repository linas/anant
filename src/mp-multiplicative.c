/**
 * mp-multiplicative.c
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

#include <mp-multiplicative.h>

// Tail-recursive helper function
static void plicplic(cpx_t result,
              void (*func)(cpx_t f, unsigned long p, int nprec),
              unsigned long  m, unsigned long n,
              int nprec)
{
	// m is prime
	if (n+1 == m)
	{
		func(result, m, nprec);
		return;
	}

	// m is not prime, and not divisible.
	if (m%n != 0)
	{
		plicplic(result, func, m, n+1, nprec);
		return;
	}
	mp_bitcnt_t bits = ((double) nprec) * 3.322 + 50;
	cpx_t fq, fn;
	cpx_init2(fq, bits);
	cpx_init2(fn, bits);

	unsigned long q = m/n;
	cpx_multiplicative(fq, func, q, nprec);
	cpx_multiplicative(fn, func, n, nprec);

	cpx_mul(result, fq, fn);
	cpx_clear(fq);
	cpx_clear(fn);
}

void cpx_multiplicative(cpx_t result,
              void (*func)(cpx_t f, unsigned long p, int nprec),
              unsigned long  n,
              int nprec)
{
	/* handle trivial case */
	if (n <= 3)
	{
		func(result, n, nprec);
		return;
	}
	plicplic(result, func, n, 2, nprec);
}

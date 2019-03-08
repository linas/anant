/**
 * mp-euler.c
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

#include "mp-binomial.h"
#include "mp-euler.h"

unsigned int cpx_euler_sum(cpx_t result,
              void (*func)(cpx_t, unsigned long, int),
              unsigned int maxterms,
              int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.322 + 50;

	mpf_t fbin;
	mpf_init2(fbin, bits);

	cpx_t term, fval;
	cpx_init2(term, bits);
	cpx_init2(fval, bits);

	mpz_t bin;
	mpz_init(bin);

	int n = 0;
	for (; n<maxterms; n++)
	{
		// `term` accumulates the inner sum
		//    $ \sum_{k=0}^n {n \choose k} f(k+1) $
		cpx_set_ui(term, 0, 0);
		for (int k=0; k<=n; k++)
		{
			func(fval, k+1, nprec);
			i_binomial_sequence(bin, n, k);
			mpf_set_z(fbin, bin);
			cpx_times_mpf(fval, fval, fbin);
			cpx_add(term, term, fval);
		}
		cpx_div_2exp(term, term, n+1);

		cpx_add(result, result, term);

		// Are we there yet?
	}

	mpz_clear(bin);
	cpx_clear(term);
	mpf_clear(fbin);
	return n;
}

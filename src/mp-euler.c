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
              unsigned int ndigits,
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

	/* Compute desired accuracy. This is used as a termination
	 * condition. We ask for twice-as-many digits, since we
	 * square the term.
	 */
	mpf_t epsi;
	mpf_init(epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(2 * 3.322*ndigits));

	mpf_t asum, aterm;
	mpf_init(asum);
	mpf_init(aterm);

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
		cpx_mod_sq(aterm, term);
		cpx_mod_sq(asum, result);
		mpf_div(aterm, aterm, asum);
		if (0 > mpf_cmp(aterm, epsi)) break;
	}

	mpf_clear(asum);
	mpf_clear(aterm);
	mpf_clear(epsi);

	mpz_clear(bin);
	cpx_clear(term);
	mpf_clear(fbin);
	return n;
}

// ============================================================
// #define TEST
#ifdef TEST

// Quick-n-dirty sunit test. Just sum $ sum_{n=0}^\infty x^n $

#include <stdio.h>

mpf_t ex;

void test_func(cpx_t f, unsigned long n, int nprec)
{
	mpf_pow_ui(f[0].re, ex, n-1);
	mpf_set_ui(f[0].im, 0);
}


int main (int argc, char * argv[])
{
   int prec, nbits;
   prec = 120;
   nbits = 3.3*prec;
   mpf_set_default_prec (nbits+200);

	double xx = 0.95;
	mpf_init2(ex, nbits);
	mpf_set_d(ex, xx);

	cpx_t result;
	cpx_init2(result, nbits);

	unsigned int nterm;
	nterm = cpx_euler_sum(result, test_func, 60, 30000000, prec);

	printf("summed to %u terms\n", nterm);

	mpf_t sum;
	mpf_init(sum);
	cpx_abs(sum, result);
	printf("got=%f expected=%f\n", mpf_get_d(sum), 1.0/(1.0-xx));
}
#endif

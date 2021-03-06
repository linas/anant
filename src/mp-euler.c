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

#include <math.h>
#include <stdbool.h>

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

	// zero out; we're accumulating into this
	cpx_set_ui(result, 0, 0);

	bool almost_done = false;
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
		// Avoid accidental loop termination by spurious zeros.
		// Look for two consecutive terms that are small.
		cpx_mod_sq(aterm, term);
		cpx_mod_sq(asum, result);
		mpf_div(aterm, aterm, asum);
		if (0 > mpf_cmp(aterm, epsi))
		{
			if (almost_done) break;
			almost_done = true;
		}
		else almost_done = false;
	}

	mpf_clear(asum);
	mpf_clear(aterm);
	mpf_clear(epsi);

	mpz_clear(bin);
	cpx_clear(term);
	cpx_clear(fval);

	mpf_clear(fbin);
	return n;
}

// ============================================================

unsigned int cpx_newton_series(cpx_t result,
              void (*func)(cpx_t, unsigned long, int),
              cpx_t zee,
              unsigned int ndigits,
              unsigned int maxterms,
              int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.322 + 50;

	mpf_t fbin;
	mpf_init2(fbin, bits);

	cpx_t term, fval, zeem1, zbin;
	cpx_init2(term, bits);
	cpx_init2(fval, bits);
	cpx_init2(zeem1, bits);
	cpx_init2(zbin, bits);

	// zeem1 = z-1
	cpx_set(zeem1, zee);
	cpx_sub_ui(zeem1, zeem1, 1, 0);

	mpz_t bin;
	mpz_init(bin);

	mpf_t asum, aterm;
	mpf_init(asum);
	mpf_init(aterm);

	/* Compute desired accuracy. This is used as a termination
	 * condition. We ask for twice-as-many digits, since we
	 * square the term.
	 */
	mpf_t epsi;
	mpf_init(epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(2 * 3.322*ndigits));

	// For integer values of zee equal to n, the binomial coefficient
	// will be exacly zero. We want to break out of the loop in such
	// cases. So...  is zee an integer, and if so, which?
	int is_int = 0;
	if (0 > mpf_cmp(zeem1[0].im, epsi))
	{
		mpf_set_d(aterm, 0.5);
		mpf_add(aterm, zeem1[0].re, aterm);
		mpf_floor(aterm, aterm);
		mpf_sub(asum, zeem1[0].re, aterm);
		mpf_mul(asum, asum, asum); // square to remove sign
		if (0 > mpf_cmp(asum, epsi))
		{
			is_int = mpf_get_ui(aterm);
		}
	}

	cpx_set_ui(result, 0, 0);

	bool almost_done = false;
	int n = 0;
	for (; n<maxterms; n++)
	{
		// If z is an integer, such that the binomial vanishes,
		// then break out of the loop.
		if (is_int && is_int == n-1) break;

		// `term` accumulates the finite difference sum
		//    $ \sum_{k=0}^n (-1)^k {n \choose k} f(k+1) $
		cpx_set_ui(term, 0, 0);
		for (int k=0; k<=n; k++)
		{
			func(fval, k+1, nprec);
			i_binomial_sequence(bin, n, k);
			mpf_set_z(fbin, bin);
			cpx_times_mpf(fval, fval, fbin);
			if (k%2 == 1) cpx_neg(fval, fval);
			cpx_add(term, term, fval);
		}

		cpx_binomial(zbin, zeem1, n);
		cpx_mul(term, term, zbin);
		if (n%2 == 1) cpx_neg(term, term);

		cpx_add(result, result, term);

		// Are we there yet?
		// Avoid accidental loop termination by spurious zeros.
		// Look for two consecutive terms that are small.
		cpx_mod_sq(aterm, term);
		cpx_mod_sq(asum, result);
		mpf_div(aterm, aterm, asum);
		if (0 > mpf_cmp(aterm, epsi))
		{
			if (almost_done) break;
			almost_done = true;
		}
		else almost_done = false;
	}

	mpf_clear(asum);
	mpf_clear(aterm);
	mpf_clear(epsi);

	mpz_clear(bin);
	cpx_clear(zeem1);
	cpx_clear(term);
	cpx_clear(fval);
	cpx_clear(zbin);
	mpf_clear(fbin);
	return n;
}

// ============================================================
// #define EULER_TEST
#ifdef EULER_TEST

// Quick-n-dirty unit test. Just sum $ sum_{n=0}^\infty x^n $

#include <stdio.h>

// Return x^n
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

// ============================================================

// #define NEWTON_TEST
#ifdef NEWTON_TEST

// Quick-n-dirty unit test. Verify interplants.

#include <stdio.h>

// Return x^n
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

	cpx_t zee, result;
	cpx_init2(zee, nbits);
	cpx_init2(result, nbits);

	mpf_t sum;
	mpf_init(sum);

	double xn = 1.0;
	for (int n=1; n<20; n++)
	{
		cpx_set_ui(zee, n, 0);
		unsigned int nterm;
		nterm = cpx_newton_series(result, test_func, zee, 20, 30000000, prec);

		printf("summed to %u terms\n", nterm);
		cpx_abs(sum, result);
		printf("n=%d got=%f expected=%f\n", n, mpf_get_d(sum), xn);

		xn *= xx;
	}
}
#endif

/*
 * mp-topsin.c
 *
 * High-precision Topologist's Sine function,
 * using the Gnu Multiple-precision library.
 *
 * Copyright (C) 2014 Linas Vepstas
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
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mp-binomial.h"
#include "mp-consts.h"
#include "mp-topsin.h"

void topsin_series (mpf_t a_k, unsigned int k, unsigned int prec)
{
	unsigned int n;
	mpf_t fourpi, numer, fact, low_bound, term, xterm;
	mpz_t bino;

	mpf_set_ui(a_k, 0);

	/* If k == 0 then we are done */
	if (0 == k) return;

	mpf_init(fourpi);
	mpf_init(numer);
	mpf_init(fact);
	mpf_init(low_bound);
	mpf_init(term);
	mpf_init(xterm);
	mpz_init(bino);

	fp_two_pi(numer, prec);
	mpf_neg(numer, numer);

	// fourpi is actually -4pi^2
	mpf_mul(fourpi, numer, numer);
	mpf_neg(fourpi, fourpi);

	mpf_set_ui(fact, 1);

	/* Get the number of binary bits from prec = log_2 10 * prec */
	long nbits = (long) floor (3.321 * prec);
	mpf_div_2exp(low_bound, fact, nbits+32);

	for (n=0; n<2023123123; n++)
	{
		i_binomial(bino, 2*n+k, 2*n);
		mpf_set_z(term, bino);
		mpf_mul(xterm, term, numer);
		mpf_div(term, xterm, fact);

		mpf_add(a_k, a_k, term);

// #define DEBUG 1
#ifdef DEBUG
		{
			double h_f, q_f, b_f;
			h_f = mpf_get_d(h);
			q_f = mpf_get_d(qmark);
			b_f = mpf_get_d(bits);
			printf("duuude place=%d, bitsdone=%ld h=%g q=%g bits=%f s=%d\n",
				place, bitsdone, h_f, q_f, b_f, mpf_sgn(h));
		}
#endif
		// If the term is small enough, we are done.
		mpf_abs(xterm, term);
		if (mpf_cmp(xterm, low_bound) < 0) break;

		// Now iterate
		mpf_mul(numer, numer, fourpi);
		mpf_mul_ui(fact, fact, (2*n+3)*2*(n+1));
	}

	if (k%2 == 0) mpf_neg(a_k, a_k);

	mpf_clear(fourpi);
	mpf_clear(numer);
	mpf_clear(fact);
	mpf_clear(low_bound);
	mpf_clear(term);
	mpf_clear(xterm);
	mpz_clear(bino);
}

/* ================================================================ */

// #define RUN_TEST
#ifdef RUN_TEST

#include "mp-trig.h"
#include <stdbool.h>

bool sum_test(double xf, int prec)
{
	bool fail = false;
	int k;
	mpf_t a_k, x, xn, term, sum, sino, low_bound;

	mpf_init(a_k);
	mpf_init(x);
	mpf_init(xn);
	mpf_init(term);
	mpf_init(sum);
	mpf_init(sino);
	mpf_init(low_bound);

	long nbits = (long) floor (3.321 * prec);
	mpf_set_ui(x, 1);
	mpf_div_2exp(low_bound, x, nbits+32);

	// Sum the series sum_k=1^\infty a_k x^k
	// It should equal sin(2pi/(1+x))
	mpf_set_d(x, xf);
	mpf_set_ui(sum, 0);
	mpf_set(xn, x);
	for (k=1; k<100000; k++)
	{
		topsin_series(a_k, k, prec);
		mpf_mul(term, a_k, xn);
		mpf_add(sum, sum, term);

		// If the term is small enough, we are done.
		mpf_abs(term, term);
		if (mpf_cmp(term, low_bound) < 0) break;

		mpf_mul(xn, xn, x);
	}

	// Now compute sin(2pi/(1+x))
	mpf_add_ui(term, x, 1);
	fp_two_pi(sino, prec);
	mpf_div(term, sino, term);

	fp_sine(sino, term, prec);

	// the sum and the sine should be equal
	mpf_sub(term, sino, sum);
	double zero = mpf_get_d(term);

	double lim = pow(10.0, -prec);
	if (fabs(zero) > lim)
	{
		printf("Error: Expecting precision 1.0e-%d got %g at x=%f\n", prec, zero, xf);
		fail = true;
	}

	mpf_clear(a_k);
	mpf_clear(x);
	mpf_clear(xn);
	mpf_clear(term);
	mpf_clear(sum);
	mpf_clear(sino);
	mpf_clear(low_bound);

	return fail;
}

int main (int argc, char * argv[])
{
	mpf_t a_k;
	int prec, nbits;

	prec = 50;

	/* Set the precision (number of binary bits) */
	/* We need more bits than what what is available, for intermediate calcs */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits+200);

	mpf_init(a_k);

	// a_1 should be -2pi
	topsin_series(a_k, 1, prec);
	double twopi = mpf_get_d(a_k);
	twopi += 2.0*M_PI;
	if (fabs(twopi) > 1.0e-16) printf("Error  at k=1: %g\n", twopi);

	// a_2 should be +2pi
	topsin_series(a_k, 2, prec);
	twopi = mpf_get_d(a_k);
	twopi -= 2.0*M_PI;
	if (fabs(twopi) > 1.0e-16) printf("Error  at k=2: %g\n", twopi);

	// a_3 should be 2pi (3-2pi^2) / 3 = -35.05851693322

	double x;
	bool fail = false;
	for (x=0.95; x>-0.95; x -= 0.018756)
	{
		bool result = sum_test(x, prec);
		fail = fail || result;
		printf("."); fflush(stdout);
	}
	printf("\n");

	if (fail) printf("Error: test failed\n");
	else printf("Success: test worked\n");

	return 0;
}
#endif

// #define PRINT_OUT_AK
#ifdef PRINT_OUT_AK
int main (int argc, char * argv[])
{
	mpf_t a_k;
	int prec, nbits;

	prec = 120;

	/* Set the precision (number of binary bits) */
	/* We need more bits than what what is available, for intermediate calcs */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits+200);

	mpf_init(a_k);

   printf("#\n# The topsin series a_k\n#\n");

	int k;
	double akprev=0.0;
	for (k=0; k<95001; k++)
	{
		topsin_series(a_k, k, prec);
		double ak = mpf_get_d(a_k);

		printf("%d	%20.16g	%20.16g\n", k, ak, ak+akprev);
		akprev = ak;
	}

	return 0;
}
#endif

/* =============================== END OF FILE =========================== */

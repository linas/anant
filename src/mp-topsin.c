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
#include "mp-topsin.h"

void topsin_series (mpf_t a_k, unsigned int k, unsigned int prec)
{
	mpf_t fourpi, numer, fact;

	mpf_set_ui(a_k, 0);

	/* If k == 0 then we are done */
	if (0 == k) return;

	mpf_init(fourpi);
	mpf_init(h);
	mpf_init(bits);
	mpf_init(one);
	mpf_init(low_bound);

	mpf_set(h, x);
	mpf_set_ui(one, 1);

	/* Get the number of binary bits from prec = log_2 10 * prec */
	long nbits = (long) floor (3.321 * prec);
	mpf_div_2exp(low_bound, one, nbits-2);

	bitsdone = -1;
	place = 1;
	while (bitsdone < nbits)
	{
		// Compute h(x) = 1/x - floor(1/x);
		mpf_ui_div(ox, 1, h);
		mpf_floor(bits, ox);
		mpf_sub(h, ox, bits);

		// bit is just floor(1/x)
		ibits = mpf_get_si(bits);
		bitsdone += ibits;

		// shift right by bits
		mpf_div_2exp(ox, one, bitsdone);

		// Add or subtract dyadic parts
		if (place%2 == 1)
		{
			mpf_add(qmark, qmark, ox);
		}
		else
		{
			mpf_sub(qmark, qmark, ox);
		}

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
		// If the remainder is zero, we are done.
		// Due to rounding precision, we can't test for explicit zero;
		// instead we test for h less than the requested precision.
		// if (0 == mpf_sgn(h)) break;
		if (mpf_cmp(h, low_bound) < 0) break;

		place ++;
	}

	mpf_clear(ox);
	mpf_clear(h);
	mpf_clear(bits);
	mpf_clear(one);
	mpf_clear(low_bound);
}

/* ================================================================ */

// #define RUN_TEST
#ifdef RUN_TEST
int main (int argc, char * argv[])
{
	int n, d;
	double qid;
	mpf_t qi, x;
	int prec, nbits;

	d = atoi(argv[1]);

	prec = 50;

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	mpf_init(qi);
	mpf_init(x);

	// n = 2*171;
	// d = 3*171;
	for (n=0; n<=d; n++)
	{
		double xd = ((double) n) / ((double) d);
		mpf_set_ui (x, n);
		mpf_div_ui (x, x, d);

		question_inverse(qi, x, prec);
		qid = mpf_get_d(qi);
		printf("%d	%f	%f\n", n, xd, qid);
	}
	
	return 0;
}
#endif

/* =============================== END OF FILE =========================== */

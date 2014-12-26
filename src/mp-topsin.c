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

	mpf_set_ui(fact, 1);

	fp_two_pi(numer, prec);
	mpf_neg(numer, numer);

	// fourpi is actually -4pi^2
	mpf_mul(fourpi, numer, numer);
	mpf_neg(fourpi, fourpi);

	/* Get the number of binary bits from prec = log_2 10 * prec */
	long nbits = (long) floor (3.321 * prec);
	mpf_div_2exp(low_bound, fact, nbits-2);

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

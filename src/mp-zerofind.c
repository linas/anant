/**
 * mp-zerofind.c
 *
 * Locate complex zeros of a function.
 *
 * Copyright (C) 2010 Linas Vepstas
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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-zerofind.h"

/* ---------------------------------------------- */
/*
 * Find bottom of parabola.  Given three points,
 * fit a parabola to that, then find the minimum.
 *
 * Return the bottom in "loc"
 */
static void quad_min(mpf_t loc,
                     mpf_t a, mpf_t b, mpf_t c,
                     mpf_t fa, mpf_t fb, mpf_t fc,
                     mp_bitcnt_t bits)
{
	mpf_t ba, bc, fba, fbc, deno, numer;
	mpf_init2 (ba, bits);
	mpf_init2 (bc, bits);
	mpf_init2 (fba, bits);
	mpf_init2 (fbc, bits);
	mpf_init2 (deno, bits);
	mpf_init2 (numer, bits);

	/* differences */
	mpf_sub (ba, b, a);
	mpf_sub (bc, b, c);
	mpf_sub (fba, fb, fa);
	mpf_sub (fbc, fb, fc);

	/* cross product */
	mpf_mul(fbc, fbc, ba);
	mpf_mul(fba, fba, bc);

	/* denominator */
	mpf_sub(deno, fbc, fba);

	/* Avoid divide by zero when denominator is exactly zero. */
	/* You'd think that this can never happen, but it does. */
	/* It would be good to throw an exception here. */
	if (mpf_sgn(deno))
	{
		/* cross again */
		mpf_mul(fbc, fbc, ba);
		mpf_mul(fba, fba, bc);

		/* numerator */
		mpf_sub(numer, fbc, fba);
		mpf_div (numer, numer, deno);
		mpf_div_ui(numer, numer, 2);

		mpf_sub(loc, b, numer);
	}
	else
	{
		mpf_set(loc, b);
	}

	mpf_clear (ba);
	mpf_clear (bc);
	mpf_clear (fba);
	mpf_clear (fbc);
	mpf_clear (deno);
	mpf_clear (numer);
}

/* =============================================== */
/**
 * cpx_find_zero.
 * Numerically locate the zero of a complex-valued function.
 *
 * @func function whose zeros are to be found.
 *       func takes z as input, returns f as output.
 *       'nprec' is the suggested decimal precision at which 'fun'
 *       should perform its calculations.
 * @initial_z initial suggestion for the location of the zero.
 * @e1, @e2 initial suggestions for a bounding ellipse. These are
 *       taken to be two vectors, specifying the major and minor axes of
 *       an ellipse, centered at 'initial_z'. The true zero is presumed
 *       to lie inside of, or at least, close to, this initial ellipse.
 * @ndigits number of decimal digits of accuracy to which the zero
 *       should be searched for.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 *
 * @returns 0 if result is valid, else an error code.
 *
 * This implements Powell's method, slightly adapted; the adaptations
 * are meant to improve convergence when 'func' is extremely noisy,
 * i.e. when any sort of quadratic behaviour is obscured by
 * high-frequency noise.
 *
 * A description of Powell's method can be found in Press, Teukolsky,
 * Vetterling, Flannery, "Numerical Recipes in C, 2nd ed.", Cambridge
 * U Press, 1999.
 *
 * This is a rather sloppy way of doing this. We really could/should
 * do something better, e.g. this:
 *    LM Delves, JN Lyness "A Numerical Method for Locating the Zeros
 *    of an Analytic Function" (1967)
 *    http://www.ams.org/journals/mcom/1967-21-100/S0025-5718-1967-0228165-4/S0025-5718-1967-0228165-4.pdf
 * or maybe this:
 *    Michael Sagraloff, Chee K. Yap, "A Simple But Exact and Efficient
 *    Algorithm for Complex Root Isolation" (2011)
 * except I'm lazy and the below mostly works.
 *
 * The below also gets used to find zeros of noisy functions: functions
 * that are not very smooth near the zero, and thus violate
 * simple-minded assumptions about analyticity.
 */

int cpx_find_zero_quad(cpx_t result,
              void (*func)(cpx_t f, cpx_t z, int nprec),
              cpx_t initial_z,
              cpx_t e1, cpx_t e2,
              int ndigits, int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.322 + 50;

	int rc = 1;
	mpf_t zero, epsi;
	mpf_init2 (zero, bits);

	/* Compute the tolerance */
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.322*ndigits));

	cpx_t s0, s1, s2, s3, sa, sb;
	cpx_init2 (s0, bits);
	cpx_init2 (s1, bits);
	cpx_init2 (s2, bits);
	cpx_init2 (s3, bits);
	cpx_init2 (sa, bits);
	cpx_init2 (sb, bits);

	cpx_t na, nb, nc;
	cpx_init2 (na, bits);
	cpx_init2 (nb, bits);
	cpx_init2 (nc, bits);

	cpx_t y0, y1, y2, y3;
	cpx_init2 (y0, bits);
	cpx_init2 (y1, bits);
	cpx_init2 (y2, bits);
	cpx_init2 (y3, bits);

	mpf_t f0, f1, f2, f3;
	mpf_init2 (f0, bits);
	mpf_init2 (f1, bits);
	mpf_init2 (f2, bits);
	mpf_init2 (f3, bits);

	mpf_t loc, lam0, lam1, lam2;
	mpf_init2 (loc, bits);
	mpf_init2 (lam0, bits);
	mpf_init2 (lam1, bits);
	mpf_init2 (lam2, bits);

	mpf_set_ui (lam0, 0);
	mpf_set_ui (lam1, 1);
	mpf_neg (lam2, lam1);

	/* Initial guess for the zero. */
	cpx_set (s0, initial_z);
	func (y0, s0, nprec);

	/* OK, so this is important, so listen-up!  The shape of the
	 * absolute value of a complex function near a zero is a 2D
	 * cone.  The sides of the cone tend to be almost straight, and
	 * so trying to fit a parabola to that ends in numerical disaster.
	 * For the quadratic minimum-finder, we really want to fit a
	 * parabola, not a cone, so that the quadratic zero finder can
	 * navigate to the bottom. Thus, squares of absolute values get
	 * used. Ugh.
	 *
	 * Of course, this is just stupid: we should take advantage of
	 * the cone-shape, and fit to that, instead.
	 */
// #define ABSVAL cpx_abs
#define ABSVAL cpx_mod_sq
	ABSVAL(f0, y0);

	/* Initial directions */
	cpx_set (na, e1);
	cpx_set (nb, e2);

	/* Iterate */
	int i;
	for (i=0; i<100; i++)
	{
		bool done1 = false;
		bool done2 = false;

		ABSVAL (zero, na);
		if (0 > mpf_cmp(zero, epsi))
		{
			cpx_set (sa, s0);
			done1 = true;
		}
		else
		{
			/* Three colinear points to start with */
			cpx_add(s1, s0, na);
			cpx_sub(s2, s0, na);

			// func (y0, s0, nprec);
			func (y1, s1, nprec);
			func (y2, s2, nprec);

			// ABSVAL(f0, y0);
			ABSVAL(f1, y1);
			ABSVAL(f2, y2);

			/* loc provides new minimum, along direction a */
			quad_min (loc, lam0, lam1, lam2, f0, f1, f2, bits);

			/* Move to that location */
			cpx_times_mpf (nc, na, loc);
			cpx_add (s3, s0, nc);
			func (y3, s3, nprec);
			ABSVAL(f3, y3);

			/* Does f3 actually improve on f0 ? */
			if (0 < mpf_cmp(f0, f3))
			{
				/* Yes, the new point is an improvement. Go there. */
				mpf_abs (zero, loc);
				if (0 < mpf_cmp(zero, epsi))
				{
					/* Not converged yet */
					cpx_times_mpf (na, na, loc);
					cpx_add (sa, s0, na);
					cpx_times_d (na, na, 0.5);
				}
				else
				{
					cpx_set(sa, s0);
					done1 = true;
				}

				/* Save last result for the next step */
				mpf_set(f0, f3);
			}
			else
			{
				bool better = false;
				/* The new point is not an improvement on f0! */
				if (0 < mpf_cmp(f0, f1))
				{
					/* ... But f1 is better */
					cpx_set (sa, s1);
					mpf_set (f0, f1);
					better = true;
				}
				if (0 < mpf_cmp(f0, f2))
				{
					/* ... But f2 is better */
					cpx_set (sa, s2);
					mpf_set (f0, f2);
					better = true;
				}
				if (better)
					cpx_times_d (na, na, 1.618);
				else
					cpx_times_d (na, na, 0.5);
			}
		}

		/* Repeat for direction b */
		ABSVAL (zero, nb);
		if (0 > mpf_cmp(zero, epsi))
		{
			cpx_set (sb, sa);
			done2 = true;
		}
		else
		{
			cpx_add(s1, sa, nb);
			cpx_sub(s2, sa, nb);

			// func (y0, sa, nprec);
			func (y1, s1, nprec);
			func (y2, s2, nprec);

			// ABSVAL(f0, y0);
			ABSVAL(f1, y1);
			ABSVAL(f2, y2);

			/* loc provides new minimum, along direction a */
			quad_min (loc, lam0, lam1, lam2, f0, f1, f2, bits);

			/* Move to that location */
			cpx_times_mpf (nc, nb, loc);
			cpx_add (s3, sa, nc);
			func (y3, s3, nprec);
			ABSVAL(f3, y3);

			/* Does f3 actually improve on f0 ? */
			if (0 < mpf_cmp(f0, f3))
			{
				/* Yes, the new point is an improvement. Go there. */
				mpf_abs (zero, loc);
				if (0 < mpf_cmp(zero, epsi))
				{
					cpx_times_mpf (nb, nb, loc);
					cpx_add (sb, sa, nb);
					cpx_times_d (nb, nb, 0.5);
				}
				else
				{
					cpx_set(sb, sa);
					done2 = true;
				}

				/* Save last result for the next step */
				mpf_set(f0, f3);
			}
			else
			{
				int better = 0;
				/* The new point is not an improvement on f0! */
				if (0 < mpf_cmp(f0, f1))
				{
					/* ... But f1 is better */
					cpx_set (sb, s1);
					mpf_set (f0, f1);
					better = 1;
				}
				if (0 < mpf_cmp(f0, f2))
				{
					/* ... But f2 is better */
					cpx_set (sb, s2);
					mpf_set (f0, f2);
					better = 1;
				}
				if (better)
					cpx_times_d (na, na, 1.618);
				else
					cpx_times_d (na, na, 0.5);
			}
		}

		/* Shuffle down */
		/* Powells' algo says shuffle, but we won't do that. */
		// cpx_set(na, nb);
		// cpx_sub(nb, sb, s0);
		cpx_set(s0, sb);

#if 0
printf("#\n# %d  done=%d %d  ", i, done1, done2);
cpx_prt("s0 = ", s0); printf("\n");
cpx_prt("# na = ", na); printf("\n");
cpx_prt("# nb = ", nb); printf("\n");
fp_prt("# min= ", f0); printf("\n");
fflush (stdout);
#endif

		if (done1 && done2)
		{
			rc = 0;
			break;
		}

#if 0
		/* Bound results away from zero.  Huh??? */
		mpf_set_d(f3, 0.05);
		if (0 > mpf_cmp(s0[0].im, f3)) break;
#endif
	}

	/* The returned value */
	cpx_set(result, s0);

	cpx_clear (s0);
	cpx_clear (s1);
	cpx_clear (s2);
	cpx_clear (s3);
	cpx_clear (sa);
	cpx_clear (sb);

	cpx_clear (na);
	cpx_clear (nb);
	cpx_clear (nc);

	cpx_clear (y0);
	cpx_clear (y1);
	cpx_clear (y2);
	cpx_clear (y3);

	mpf_clear (f0);
	mpf_clear (f1);
	mpf_clear (f2);
	mpf_clear (f3);

	mpf_clear (loc);
	mpf_clear (lam0);
	mpf_clear (lam1);
	mpf_clear (lam2);

	mpf_clear (zero);
	mpf_clear (epsi);

	return rc;
}

/* =============================================== */
// Try to fit a conic section.
//
// Assume that the the points fa, fb, fc sit on the surface of a cone,
// and that the tip of the cone is the zero. Just fit the conic section
// and use that to estimate the location of the tip of the cone.

static void conic(cpx_t loc,
                  cpx_t za, cpx_t zb, cpx_t zc,
                  cpx_t fa, cpx_t fb, cpx_t fc,
                  mp_bitcnt_t bits)
{
	cpx_t df;
	cpx_init2(df, bits);

	cpx_t zba, zca;
	cpx_init2(zba, bits);
	cpx_init2(zca, bits);

	cpx_sub(zba, zb, za);
	cpx_sub(zca, zc, za);

	mpf_t lb, lc;
	mpf_init2(lb, bits);
	mpf_init2(lc, bits);

	cpx_abs(lb, zba);
	cpx_abs(lc, zca);

	// This is true, if lb is longer than lc
	// Always use the longer arm of extrapolation.
	if (0 < mpf_cmp(lb, lc))
	{
		cpx_sub(df, fb, fa);
		cpx_div(loc, zba, df);
	}
	else
	{
		cpx_sub(df, fc, fa);
		cpx_div(loc, zca, df);
	}

	cpx_mul(loc, loc, fa);
	cpx_sub(loc, za, loc);

	mpf_clear(lb);
	mpf_clear(lc);
	cpx_clear(zba);
	cpx_clear(zca);
	cpx_clear(df);
}

/* =============================================== */
/**
 * cpx_find_zero.
 * Numerically locate the zero of a complex-valued function.
 *
 * @func function whose zeros are to be found.
 *       func takes z as input, returns f as output.
 *       'nprec' is the suggested decimal precision at which 'fun'
 *       should perform its calculations.
 * @initial_z initial suggestion for the location of the zero.
 * @e1, @e2 initial suggestions for a bounding ellipse. These are
 *       taken to be two vectors, specifying the major and minor axes of
 *       an ellipse, centered at 'initial_z'. The true zero is presumed
 *       to lie inside of, or at least, close to, this initial ellipse.
 * @ndigits number of decimal digits of accuracy to which the zero
 *       should be searched for.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 *
 * @returns 0 if result is valid, else an error code.
 *
 * This implements a brute-force conic-section fit.
 *
 * This is a rather sloppy way of doing this. We really could/should
 * do something better, e.g. this:
 *    LM Delves, JN Lyness "A Numerical Method for Locating the Zeros
 *    of an Analytic Function" (1967)
 *    http://www.ams.org/journals/mcom/1967-21-100/S0025-5718-1967-0228165-4/S0025-5718-1967-0228165-4.pdf
 * or maybe this:
 *    Michael Sagraloff, Chee K. Yap, "A Simple But Exact and Efficient
 *    Algorithm for Complex Root Isolation" (2011)
 * except I'm lazy and the below mostly works.
 */

int cpx_find_zero(cpx_t result,
              void (*func)(cpx_t f, cpx_t z, int nprec),
              cpx_t initial_z,
              cpx_t e1, cpx_t e2,
              int ndigits, int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.322 + 50;

	int rc = 1;
	mpf_t zero, epsi;
	mpf_init2 (zero, bits);

	/* Compute the tolerance */
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.322*ndigits));

	cpx_t s0, s1, s2, s3;
	cpx_init2 (s0, bits);
	cpx_init2 (s1, bits);
	cpx_init2 (s2, bits);
	cpx_init2 (s3, bits);

	cpx_t y0, y1, y2, y3;
	cpx_init2 (y0, bits);
	cpx_init2 (y1, bits);
	cpx_init2 (y2, bits);
	cpx_init2 (y3, bits);

	mpf_t f0, f1, f2, f3;
	mpf_init2 (f0, bits);
	mpf_init2 (f1, bits);
	mpf_init2 (f2, bits);
	mpf_init2 (f3, bits);

	/* Initial guess for the zero. */
	cpx_set (s0, initial_z);
	cpx_add (s1, initial_z, e1);
	cpx_add (s2, initial_z, e2);

	func (y0, s0, nprec);
	func (y1, s1, nprec);
	func (y2, s2, nprec);

	cpx_abs(f0, y0);
	cpx_abs(f1, y1);
	cpx_abs(f2, y2);

	// Place into sorted order.
	if (0 < mpf_cmp(f0, f1))
	{
		cpx_set(y3, y0); cpx_set(y0, y1); cpx_set(y1, y3);
		cpx_set(s3, s0); cpx_set(s0, s1); cpx_set(s1, s3);
		mpf_set(f3, f0); mpf_set(f0, f1); mpf_set(f1, f3);
	}
	if (0 < mpf_cmp(f0, f2))
	{
		cpx_set(y3, y0); cpx_set(y0, y2); cpx_set(y2, y3);
		cpx_set(s3, s0); cpx_set(s0, s2); cpx_set(s2, s3);
		mpf_set(f3, f0); mpf_set(f0, f2); mpf_set(f2, f3);
	}
	if (0 < mpf_cmp(f1, f2))
	{
		cpx_set(y3, y1); cpx_set(y1, y2); cpx_set(y2, y3);
		cpx_set(s3, s1); cpx_set(s1, s2); cpx_set(s2, s3);
		mpf_set(f3, f1); mpf_set(f1, f2); mpf_set(f2, f3);
	}

	/* Iterate */
	int i;
	for (i=0; i<100; i++)
	{
		cpx_sub(s3, s1, s0);
		cpx_abs (zero, s3);
		if (0 > mpf_cmp(zero, epsi))
		{
			rc = 0;
			cpx_set (result, s0);
			break;
		}

		conic(s3, s0, s1, s2, y0, y1, y2, bits);
		func (y3, s3, nprec);
		cpx_abs(f3, y3);

		if (0 < mpf_cmp(f3, f2))
		{
			// Oh no Mr bill!  Its not an improvement!
			cpx_set (result, s0);
			rc = 1;
			break;
		}
		else
		{
			// Copy over
			cpx_set(s2, s3);
			cpx_set(y2, y3);
			mpf_set(f2, f3);

			// Place into sorted order.
			if (0 < mpf_cmp(f0, f1))
			{
				cpx_set(y3, y0); cpx_set(y0, y1); cpx_set(y1, y3);
				cpx_set(s3, s0); cpx_set(s0, s1); cpx_set(s1, s3);
				mpf_set(f3, f0); mpf_set(f0, f1); mpf_set(f1, f3);
			}
			if (0 < mpf_cmp(f0, f2))
			{
				cpx_set(y3, y0); cpx_set(y0, y2); cpx_set(y2, y3);
				cpx_set(s3, s0); cpx_set(s0, s2); cpx_set(s2, s3);
				mpf_set(f3, f0); mpf_set(f0, f2); mpf_set(f2, f3);
			}
			if (0 < mpf_cmp(f1, f2))
			{
				cpx_set(y3, y1); cpx_set(y1, y2); cpx_set(y2, y3);
				cpx_set(s3, s1); cpx_set(s1, s2); cpx_set(s2, s3);
				mpf_set(f3, f1); mpf_set(f1, f2); mpf_set(f2, f3);
			}
		}

#if 0
cpx_sub(s3, s1, s0);
cpx_abs (zero, s3);
printf("#\n# %d ", i);
cpx_prt("s0 = ", s0); printf("\n");
fp_prt("# err= ", zero); printf("\n");
fp_prt("# min= ", f0); printf("\n");
printf("\n");
fflush (stdout);
#endif

	}

	cpx_clear (s0);
	cpx_clear (s1);
	cpx_clear (s2);
	cpx_clear (s3);

	cpx_clear (y0);
	cpx_clear (y1);
	cpx_clear (y2);
	cpx_clear (y3);

	mpf_clear (f0);
	mpf_clear (f1);
	mpf_clear (f2);
	mpf_clear (f3);

	mpf_clear (zero);
	mpf_clear (epsi);

	return rc;
}

/* =============================================== */

// #define TEST
#ifdef TEST
int iter = 0;
static void foo(cpx_t y, cpx_t s, int nprec)
{
cpx_prt("called with= ", s); printf("\n");
	cpx_t cent;
	cpx_init (cent);
	cpx_set_d(cent, 0.4980812345, 18.313412345);
	cpx_sub(cent, s, cent);
	cpx_times_d(y, cent, 5.11);
	cpx_mul(cent, cent, cent);
	cpx_times_d(cent, cent, 35.11);
	cpx_add(y, y, cent);
cpx_prt("returning ", y); printf("\n");
	cpx_clear (cent);
	iter ++;
}

int main()
{
	cpx_t z0, zg, e1, e2;
	cpx_init (z0);
	cpx_init (zg);
	cpx_init (e1);
	cpx_init (e2);

	cpx_set_d(zg, 0.5, 18.0);
	cpx_set_d(e1, 1.0, 0.0);
	cpx_set_d(e2, 0.0, 1.0);

	iter = 0;
	int rc = cpx_find_zero(z0, foo, zg, e1, e2, 16, 20);

	printf ("conic got rc=%d after %d iter\n", rc, iter);
	cpx_prt("z0 = ", z0);
	printf("\n\n\n\n");

#if 1
	iter = 0;
	rc = cpx_find_zero_quad(z0, foo, zg, e1, e2, 16, 20);

	printf ("quad got rc=%d after %d iter\n", rc, iter);
	cpx_prt("z0 = ", z0);
#endif

	printf("\n");
}

#endif

/* =============================================== */

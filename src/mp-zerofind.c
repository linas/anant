/**
 * mp-zerofind.c
 *
 * Locate zeros of a function. 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mp-complex.h"
#include "mp-zerofind.h"

#ifdef TEST
static void test_parabola(cpx_t y, cpx_t s, int nprec)
{
	cpx_t cent;
	cpx_init (cent);
	cpx_set_d(cent, 0.4980812345, 18.313412345);
	cpx_sub(cent, cent, s);
	cpx_mul(y, cent, cent);
	cpx_clear (cent);
}
#endif

/* ---------------------------------------------- */
/*
 * Find bottom of parabola.  Given three points,
 * fit a parabola to that, then find the minimum.
 * 
 * return the bottom in "loc"
 */
static void quad_min(mpf_t loc, mpf_t a, mpf_t b, mpf_t c,
              mpf_t fa, mpf_t fb, mpf_t fc)
{
	mpf_t ba, bc, fba, fbc, deno, numer;
	mpf_init (ba);
	mpf_init (bc);
	mpf_init (fba);
	mpf_init (fbc);
	mpf_init (deno);
	mpf_init (numer);

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

	/* cross again */
	mpf_mul(fbc, fbc, ba);
	mpf_mul(fba, fba, bc);

	/* numerator */
	mpf_sub(numer, fbc, fba);
	mpf_div (numer, numer, deno);
	mpf_div_ui(numer, numer, 2);

	mpf_sub(loc, b, numer);

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
 *       'nprec' is the suggested decimal precisiton at which 'fun'
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
 * are meant to imporve convergence when 'func' is extremely noisy, 
 * i.e. when any sort of quardatic behaviour is obscured by
 * high-freqency noise.
 *
 * A description of Powell's method can be found in Press, Teukolsky,
 * Vetterling, Flannery, "Numerical Recipes in C, 2nd ed.", Cambridge 
 * U Press, 1999.
 */

int cpx_find_zero(cpx_t result,
              void (*func)(cpx_t f, cpx_t z, int nprec),
              cpx_t initial_z,
              cpx_t e1, cpx_t e2,
              int ndigits, int nprec)
{
	int rc = 1;
	mpf_t zero, epsi;
	mpf_init (zero);

	/* compute the tolerance */
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*ndigits));

	cpx_t s0, s1, s2, s3, sa, sb;
	cpx_init (s0);
	cpx_init (s1);
	cpx_init (s2);
	cpx_init (s3);
	cpx_init (sa);
	cpx_init (sb);

	cpx_t na, nb, nc;
	cpx_init (na);
	cpx_init (nb);
	cpx_init (nc);

	cpx_t y0, y1, y2, y3;
	cpx_init (y0);
	cpx_init (y1);
	cpx_init (y2);
	cpx_init (y3);

	mpf_t f0, f1, f2, f3;
	mpf_init (f0);
	mpf_init (f1);
	mpf_init (f2);
	mpf_init (f3);

	mpf_t loc, lam0, lam1, lam2;
	mpf_init (loc);
	mpf_init (lam0);
	mpf_init (lam1);
	mpf_init (lam2);

	mpf_set_ui (lam0, 0);
	mpf_set_ui (lam1, 1);
	mpf_neg (lam2, lam1);

	/* Initial guess for the zero. */
	cpx_set (s0, initial_z);
	func (y0, s0, nprec);
	cpx_abs(f0, y0);

	/* Initial directions */
	cpx_set (na, e1);
	cpx_set (nb, e2);

	/* Iterate */
	int i;
	for (i=0; i<100; i++)
	{
		int done1 = 0;
		int done2 = 0;

		cpx_abs (zero, na);
		if (0 > mpf_cmp(zero, epsi))
		{
			cpx_set (sa, s0);
			done1 = 1;
		}
		else
		{
			/* Three colinear points to start with */
			cpx_add(s1, s0, na);
			cpx_sub(s2, s0, na);

			// func (y0, s0, nprec);
			func (y1, s1, nprec);
			func (y2, s2, nprec);

			// cpx_abs(f0, y0);
			cpx_abs(f1, y1);
			cpx_abs(f2, y2);

			/* loc provides new minimum, along direction a */
			quad_min (loc, lam0, lam1, lam2, f0, f1, f2);

			/* Move to that location */
			cpx_times_mpf (nc, na, loc);
			cpx_add (s3, s0, nc);
			func (y3, s3, nprec);
			cpx_abs(f3, y3);

			/* Does f3 actually improve on f0 ? */
			if (0 < mpf_cmp(f0, f3))
			{
				/* Yes, the new point is an improvement. Go there. */
				mpf_abs (zero, loc);
				if (0 < mpf_cmp(zero, epsi))
				{
					cpx_times_mpf (na, na, loc);
					cpx_add (sa, s0, na);
					cpx_times_d (na, na, 0.5);
				}
				else
				{
					cpx_set(sa, s0);
					done1 = 1;
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
					cpx_set (sa, s1);
					mpf_set (f0, f1);
					better = 1;
				}
				if (0 < mpf_cmp(f0, f2))
				{
					/* ... But f2 is better */
					cpx_set (sa, s2);
					mpf_set (f0, f2);
					better = 1;
				}
				if (better)
					cpx_times_d (na, na, 1.618);
				else
					cpx_times_d (na, na, 0.5);
			}
		}

		/* Repeat for direction b */
		cpx_abs (zero, nb);
		if (0 > mpf_cmp(zero, epsi))
		{
			cpx_set (sb, sa);
			done2 = 1;
		}
		else
		{
			cpx_add(s1, sa, nb);
			cpx_sub(s2, sa, nb);

			// func (y0, sa, nprec);
			func (y1, s1, nprec);
			func (y2, s2, nprec);

			// cpx_abs(f0, y0);
			cpx_abs(f1, y1);
			cpx_abs(f2, y2);

			/* loc provides new minimum, along direction a */
			quad_min (loc, lam0, lam1, lam2, f0, f1, f2);

			/* Move to that location */
			cpx_times_mpf (nc, nb, loc);
			cpx_add (s3, sa, nc);
			func (y3, s3, nprec);
			cpx_abs(f3, y3);

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
					done2 = 1;
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
printf("#\n# %d  ", i);
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

		/* bound results away from zero */
		mpf_set_ui(f3, 1);
		mpf_div_ui(f3, f3, 20);
		if (0 > mpf_cmp(s0[0].im, f3)) break;	
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


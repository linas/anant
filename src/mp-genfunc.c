/*
 * Generating functions for assorted number-theoretic functions
 * Implementation in bignums.
 *
 * Linas Vepstas - April 2016, October 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-trig.h>
#include <mp-complex.h>
#include <mp-consts.h>

#include "genfunc.h"

/*
 * Ordinary generating function for function func.
 * Computes ogf(z) = sum_{n=1}^\infty func(n) z^n
 *
 * Assumes that func(n) is some arithmetic series.
 * Assumes that ogf(z) does not converge for |z|>=1
 * Assumes that ogf9z) converges poorly near |z|=1, so that
 *    some termination measures are taken, so that the sum does
 *    not run forever. The termination measures assume that
 *    func(n) is bounded by n.  i.e. this is to gaurantee good
 *    data when the system is not overflowing.
 */
void cpx_ordinary_genfunc(cpx_t sum, cpx_t z, int prec, long (*func)(long))
{
	mpf_t zabs, gabs, epsi;
	mpf_init (gabs);
	mpf_init (zabs);
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	// Not defined for |z| > 1
	mpf_sub_ui(gabs, gabs, 1);
	mpf_neg(gabs, gabs);
	if (0 > mpf_cmp(gabs, epsi)) return;

	// Limit the number of iterations as we approach the edge.
	// This assumes that func(n) is bounded by n
	double dist_to_circle = mpf_get_d(gabs);
	int niter = ceil (2.302585*prec / dist_to_circle);
	niter += ceil (log(niter) / dist_to_circle); // assume func bounded by n

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	for (int n=1; n < niter ; n++)
	{
		long funv = func(n);
		if (0 != funv)
		{
			cpx_times_ui(term, zn, labs(funv));
			if (funv < 0)
				cpx_sub(sum, sum, term);
			else
				cpx_add(sum, sum, term);
		}
		cpx_mul(zn, zn, z);

#if SLOW_VERSION_NOT_USING_NITER
		// The following checks the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, zn);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) return;
#endif
	}
}


/**
 * Exponential generating function for the arithmetic series
 *
 * Computes egf(z) = exp(-|z|) sum_{n=1}^infty func(n) z^n / n!
 *
 * Note the assumption about the leading asymptotic behavior of
 * the series.
 */
void cpx_exponential_genfunc(cpx_t sum, cpx_t z, int prec, long (*func)(long))
{
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;

	mpf_t zabs, gabs, epsi, fact;
	mpf_init2 (gabs, bits);
	mpf_init2 (zabs, bits);
	mpf_init2 (epsi, bits);
	mpf_init2 (fact, bits);
	mpf_set_ui(fact, 1);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.322*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	cpx_t zn, term;
	cpx_init2(zn, bits);
	cpx_init2(term, bits);
	cpx_set(zn, z);

	for (int n=1; ; n++)
	{
		long funv = func(n);
		if (0 != funv)
		{
			cpx_times_ui(term, zn, labs(funv));
			cpx_times_mpf(term, term, fact);
			if (funv < 0)
				cpx_sub(sum, sum, term);
			else
				cpx_add(sum, sum, term);

			// The following checks the loop termination condition,
			// which is that the size of the term is less than epsilon.
			cpx_abs(gabs, term);
			mpf_mul_ui(gabs, gabs, n);

			cpx_abs(zabs, sum);
			mpf_mul(zabs, zabs, epsi);

			// if (n * zn/n! < epsi * sum) return;
			if (0 > mpf_cmp(gabs, zabs)) break;
		}

		cpx_mul(zn, zn, z);
		mpf_div_ui(fact, fact, n+1);
	}

	// Remove the leading exponential order.
	cpx_abs(gabs, z);
	mpf_neg(gabs, gabs);
	fp_exp(gabs, gabs, prec);

	cpx_times_mpf(sum, sum, gabs);
}

/**
 * Exponential generating function for the arithmetic series.
 * Same as above, except that func returns an real value, stored
 * in mpf_t
 */
void cpx_exponential_genfunc_mpf(cpx_t sum, cpx_t z, int prec,
                                 void (*func)(mpf_t, long))
{
	mpf_t zabs, gabs, epsi, fact, func_val, fabs;
	mpf_init (gabs);
	mpf_init (zabs);
	mpf_init (epsi);
	mpf_init (fact);
	mpf_init (func_val);
	mpf_init (fabs);
	mpf_set_ui(fact, 1);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	int n;
	for (n=1; ; n++)
	{
		func(func_val, n);

		// The below is a weird hack to check for a value of zero
		// returned by func.  It fails, if func is trying to return
		// a very small but non-zero value; that is why it a hack.
		// Currently, none of our functions return small values,
		// so this is OK, for now.
		mpf_abs(fabs, func_val);
		if (0 < mpf_cmp(fabs, epsi))
		{
			cpx_times_mpf(term, zn, func_val);
			cpx_times_mpf(term, term, fact);
			cpx_add(sum, sum, term);

			// The following checks the loop termination condition,
			// which is that the size of the term is less than epsilon.
			cpx_abs(gabs, term);
			// mpf_mul_ui(gabs, gabs, n);
			mpf_mul(gabs, gabs, fabs);

			cpx_abs(zabs, sum);
			mpf_mul(zabs, zabs, epsi);

			// if (n * zn/n! < epsi * sum) return;
			if (0 > mpf_cmp(gabs, zabs)) break;
		}

		cpx_mul(zn, zn, z);
		mpf_div_ui(fact, fact, n+1);
	}

	// Remove the leading exponential order.
	cpx_abs(gabs, z);
	mpf_neg(gabs, gabs);
	fp_exp(gabs, gabs, prec);

	cpx_times_mpf(sum, sum, gabs);

#if 0
double s = mpf_get_d(gabs);
cpx_abs(gabs, z);
double r = mpf_get_d(gabs);
cpx_abs(gabs, sum);
double g = mpf_get_d(gabs);
double ph = 0.5 + 0.5 * atan2(cpx_get_im(z), cpx_get_re(z))/M_PI;
printf("duuude r=%9.3f ph=%f n=%d g=%9.5e scale=%g\n", r, ph, n, g,s);
#endif
	mpf_clear (gabs);
	mpf_clear (zabs);
	mpf_clear (epsi);
	mpf_clear (fact);
	mpf_clear (func_val);
	mpf_clear (fabs);
}

/**
 * Twist the z value by gamma.
 */
void cpx_exponential_twist(cpx_t sum, cpx_t z, int prec, long (*func)(long))
{
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;

	cpx_t zt;
	cpx_init2(zt, bits);

	mpf_t r, t, g, gd, pi;
	mpf_init2(r, bits);
	mpf_init2(t, bits);
	mpf_init2(g, bits);
	mpf_init2(gd, bits);
	mpf_init2(pi, bits);
	fp_pi(pi, prec);

	// Decompose z = r exp itheta
	cpx_abs(r, z);
	fp_arctan2(t, z->im, z->re, prec);
	mpf_div(t, t, pi);

	// t runs between -1 and 1.
	int sgn = mpf_sgn (t);
	if (sgn < 0)
	{
		mpf_neg(t, t);
	}

	// rescale so that t runs between 0 and 1/2
	mpf_div_ui(t, t, 2);

	// ---------------------------------------
	// Compute gamma = g(t) = t/(n+t)

	unsigned long n = 1;
	mpf_add_ui (g, t, n);
	mpf_div (g, t, g);

	// ---------------------------------------
	// rescale r as well. Need g-dot(t) = n/(n+t)^2
	mpf_add_ui (gd, t, n);
	mpf_ui_div (gd, 1, gd);
	mpf_mul (gd, gd, gd);
	mpf_mul_ui (gd, gd, n);

	// --------------------------------------
	// g runs between 0 and 1/3 = 1/(n+2)
	mpf_set(t, g);

	// rescale back up that t is doubled again
	mpf_mul_ui(t, t, 2);

	// Multiply by pi, so that t runs over a sub-range of
	// -pi to pi.  For a single power of g, that means that
	// t runs over 0 to 2pi/3 = 2pi/(n+2)
	mpf_mul(t, t, pi);

	// --------------------------------------
	// Rescale -- take power r^gdot
	fp_log (r, r, prec);
	mpf_mul (r, r, gd);
	fp_exp(r, r, prec);

	// --------------------------------------
	// reconstruct z.
	// zt = r exp (it)
	cpx_set_ui(zt, 0, 1);           // zt = i
	cpx_times_mpf(zt, zt, t);       // zt = itheta
	cpx_exp(zt, zt, prec);          // zt = exp(itheta)
	cpx_times_mpf(zt, zt, r);       // zt = r exp(itheta)

	cpx_exponential_genfunc(sum, zt, prec, func);

	cpx_clear(zt);
	mpf_clear(r);
	mpf_clear(t);
	mpf_clear(g);
	mpf_clear(gd);
	mpf_clear(pi);
}

/*
 * mp-trig.c
 *
 * High-precison Elementary functions, using the 
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2005, 2006, 2010 Linas Vepstas
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
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* ======================================================================= */
/**
 * i_pow - raise n to the m power
 */

void i_pow (mpz_t p, unsigned int n, unsigned int m)
{
	DECLARE_I_CACHE (cache);
	cache.disabled = 1;

	if ((1 == n) || (0 == m))
	{
		mpz_set_ui (p, 1); 
		return;
	}

	int hit = i_triangle_cache_check (&cache, n+m, m);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, p, n+m, m);
	}
	else
	{
		i_pow (p, n, m-1);
		mpz_mul_ui (p, p, n);
		i_triangle_cache_store (&cache, p, n+m, m);
	}
}

/**
 * fp_inv_pow - raise n to the -m power, where m must be positive. 
 */

void fp_inv_pow (mpf_t p, unsigned int n, unsigned int m)
{
	DECLARE_FP_CACHE (cache);
	if (1 == n)
	{
		mpf_set_ui (p, 1); 
		return;
	}

	int hit = fp_triangle_cache_check (&cache, n+m, m);
	if (hit)
	{
		fp_triangle_cache_fetch (&cache, p, n+m, m);
	}
	else
	{
		mpz_t ip;
		mpz_init (ip);
		i_pow (ip, n, m);
		mpf_set_z (p, ip);
		mpf_ui_div (p, 1, p);
		mpz_clear (ip);
		fp_triangle_cache_store (&cache, p, n+m, m, 1);
	}
}

/* ======================================================================= */
/**
 * fp_exp -  Floating point exponential
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

/* exp_helpter is not static, because its sneakily used by fp_e to return e */
void fp_exp_helper (mpf_t ex, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_set (z_n, zee);
	
	mpf_set_ui (ex, 1);
	mpf_set_ui (fact, 1);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);

	unsigned int n=1;
	while(1)
	{
		fp_inv_factorial (fact, n, prec);
		mpf_mul (term, z_n, fact);
		mpf_add (ex, ex, term);
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, zee);
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (maxterm);
}

/**
 * fp_exp -- return exp(z) for real z
 *
 * Uses two-part algo for computing this:
 * It computes z = n+x so that |x|<0.5
 * Thus, e^z = e^n e^x
 * e^n is quickly computed as an integer power.
 * e^x is computed by classical taylor's series.
 */
void fp_exp (mpf_t ex, const mpf_t z, unsigned int prec)
{
	mpf_t zee, zf;
	mpf_init (zee);
	mpf_init (zf);
	
	mpf_set (zee, z);
	mpf_floor (zf, zee);
	mpf_sub (zee, zee, zf);
	long intpart = mpf_get_si (zf);

	if (mpf_cmp_d (zee, 0.5) > 0)
	{
		mpf_sub_ui (zee, zee, 1);
		intpart ++;
	}
	fp_exp_helper (ex, zee, prec);

	/* If there's an integerpart, compute e^intpart */
	if (intpart)
	{
		fp_e (zee, prec);
		if (0 < intpart)
		{
			mpf_pow_ui (zf, zee, intpart);
			mpf_mul (ex, ex, zf);
		}
		else
		{
			mpf_pow_ui (zf, zee, -intpart);
			mpf_div (ex, ex, zf);
		}
	}

	mpf_clear (zee);
	mpf_clear (zf);
}

/* ======================================================================= */
/**
 * fp_sine_series -  Floating point sine function
 * Implemented using a brute-force, very simple algo: 
 * i.e. sum the Taylor expansion. No other adjustments made.
 */

static void fp_sine_series (mpf_t si, const mpf_t z, unsigned int prec)
{
	mpf_t zsq, z_n, fact, term;

	mpf_init (zsq);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (z_n, z);
	mpf_mul (zsq, z, z);
	mpf_set_ui (si, 0);
	
	/* Use 10^{-prec} for smallest term in sum */
	prec += 2;
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);

	unsigned int n=1;
	unsigned int s=0;
	while(1)
	{
		fp_inv_factorial (fact, n, prec);
		mpf_mul (term, z_n, fact);

		if (0 == s%2)
		{
			mpf_add (si, si, term);
		}
		else
		{
			mpf_sub (si, si, term);
		}
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n+=2;
		mpf_mul (z_n, z_n, zsq);

		s++;
	}
	
	mpf_clear (zsq);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (maxterm);
}

/**
 * fp_sine -  Floating point sine function
 *
 * Adds or subtracts multiples of pi-halves until the argument
 * is in the first or fourth quadrant, and then performs the
 * brute-force taylor summation on it. 
 */

void fp_sine (mpf_t si, const mpf_t z, unsigned int prec)
{
	mpf_t zee, pih, per, top;

	mpf_init (zee);
	mpf_init (pih);
	mpf_init (per);
	mpf_init (top);

	/* Make copy of argument now! */
	mpf_set (zee, z);

	/* subtract off multiple of pi-halves */
	fp_two_over_pi (top, prec);
	mpf_mul (per, zee, top);
	mpf_floor (per, per);
	long quad = mpf_get_si (per);
	
	fp_pi_half (pih, prec);
	mpf_mul (per, per, pih);
	mpf_sub (zee, zee, per);

	/* adjust for the quadrant */
	unsigned long iq = labs (quad);
	if ((1 == iq%4) || (3 == iq%4))
	{
		mpf_sub (zee, pih, zee);
	}
	fp_sine_series (si, zee, prec);

	/* bump for negative quadrants */
	if (0>quad) iq++;
	iq /= 2;
	if (iq%2)
	{
		mpf_neg (si, si);
	}

	mpf_clear (zee);
	mpf_clear (pih);
	mpf_clear (per);
	mpf_clear (top);
}

/* ======================================================================= */
#if OBSOLETE_FOR_REFERENCE_ONLY
/* The cosine_series code below works great, except that it has poorer
 * convergence for large z values. Thus, its re-implemented in terms 
 * of sine, which is more efficient.
 *
 * This routine should be useful as an alternate high-precision algo, 
 * useful for testing/validation, and should be moved to the test package.
 */

static void fp_cosine_series (mpf_t co, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_mul (z_n, zee, zee);
	mpf_set_ui (co, 1);
	mpf_set_ui (fact, 2);
	
	/* Use 10^{-prec} for smallest term in sum */
	prec +=2;
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);

	unsigned int n=2;
	unsigned int s=1;
	while(1)
	{
		mpf_div (term, z_n, fact);

		if (0 == s%2)
		{
			mpf_add (co, co, term);
		}
		else
		{
			mpf_sub (co, co, term);
		}
		
		/* Don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);

		s++;
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (maxterm);
}
#endif /* OBSOLETE_FOR_REFERENCE_ONLY */

/**
 * fp_cosine -  Floating point cosine function
 * Implemented on top of the sine function, just to keep things simple, stupid.
 */
void fp_cosine (mpf_t si, const mpf_t z, unsigned int prec)
{
	mpf_t zee, pih;

	mpf_init (zee);
	mpf_init (pih);

	/* Make copy of argument now! */
	mpf_set (zee, z);

	/* add pi-halves */
	fp_pi_half (pih, prec);
	mpf_add (zee, zee, pih);
	fp_sine(si, zee, prec);

	mpf_clear (zee);
	mpf_clear (pih);
}

/* ======================================================================= */
/**
 * cpx_exp -  complex-valued exp, built out of the real-valued funcs.
 */

void cpx_exp (cpx_t ex, const cpx_t z, unsigned int prec)
{
	mpf_t mag, si, co;

	mpf_init (mag);
	mpf_init (si);
	mpf_init (co);

	fp_exp (mag, z->re, prec);
	fp_cosine (co, z->im, prec);
	fp_sine (si, z->im, prec);

	mpf_mul (ex->re, mag, co);
	mpf_mul (ex->im, mag, si);
	
	mpf_clear (mag);
	mpf_clear (si);
	mpf_clear (co);
}

void cpx_sine (cpx_t sn, const cpx_t z, unsigned int prec)
{
	cpx_t zee;
	cpx_init (zee);
	cpx_times_i (zee, z);
	
	cpx_exp (sn, zee, prec);
	cpx_neg (zee, zee);
	cpx_exp (zee, zee, prec);
	cpx_sub (sn, sn, zee);

	/* divide by two; using mpf_div_2exp() is faster than mpf_div_ui() */
	// cpx_div_ui (sn, sn, 2);
	mpf_div_2exp (sn[0].re, sn[0].re, 1);
	mpf_div_2exp (sn[0].im, sn[0].im, 1);
	cpx_times_i (sn, sn);
	cpx_neg (sn, sn);
	
	cpx_clear (zee);
}

void cpx_cosine (cpx_t cs, const cpx_t z, unsigned int prec)
{
	cpx_t zee;
	cpx_init (zee);
	cpx_times_i (zee, z);
	
	cpx_exp (cs, zee, prec);
	cpx_neg (zee, zee);
	cpx_exp (zee, zee, prec);
	cpx_add (cs, cs, zee);
	cpx_div_ui (cs, cs, 2);
	
	cpx_clear (zee);
}

void cpx_tangent (cpx_t tn, const cpx_t z, unsigned int prec)
{
	cpx_t zee, cn, sn;
	cpx_init (zee);
	cpx_init (cn);
	cpx_init (sn);
	
	cpx_times_i (zee, z);
	cpx_exp (tn, zee, prec);
	cpx_neg (zee, zee);
	cpx_exp (sn, zee, prec);
	
	/* cn now holds 2*cos */
	cpx_add (cn, tn, sn);
	
	/* sn now holds 2i*sin */
	cpx_sub (sn, tn, sn);
	
	/* tan = sin/cos */
	cpx_div (tn, sn, cn);
	cpx_times_i (tn, tn);
	
	cpx_clear (zee);
	cpx_clear (sn);
	cpx_clear (cn);
}

/* ======================================================================= */
/**
 * fp_log_m1 -  Floating point logarithm
 * Computes -log(1-z) using Taylor expansion for small z.
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_log_m1 (mpf_t lg, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_mul (z_n, zee, zee);
	mpf_set (lg, zee);
	
	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);

	unsigned int n=2;
	while(1)
	{
		mpf_div_ui (term, z_n, n);
		mpf_add (lg, lg, term);
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n ++;
		mpf_mul (z_n, z_n, zee);
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (term);

	mpf_clear (maxterm);
}

static inline void fp_log_simple (mpf_t lg, const mpf_t z, unsigned int prec)
{
	mpf_t zee;
	mpf_init (zee);
	if (mpf_cmp_d(z, 1.618) > 0)
	{
		mpf_ui_div (zee, 1, z);
		mpf_ui_sub (zee, 1, zee);
		fp_log_m1 (lg, zee, prec);
	}
	else
	{
		mpf_ui_sub (zee, 1, z);
		fp_log_m1 (lg, zee, prec);
		mpf_neg (lg, lg);
	}
	mpf_clear (zee);
}

/* compute and cache logarithm of 1+2^-k */
void fp_log_2exp (mpf_t lg, unsigned int k, unsigned int prec)
{
	DECLARE_FP_CACHE (log_n);

	if (prec <= fp_one_d_cache_check (&log_n, k))
	{
		fp_one_d_cache_fetch (&log_n, lg, k);
		return;
	}
	
	mpf_set_ui (lg, 1);
	mpf_div_2exp (lg, lg, k);
	mpf_neg (lg, lg);
	fp_log_m1 (lg, lg, prec);
	mpf_neg (lg, lg);

	fp_one_d_cache_store (&log_n, lg, k, prec);
}

/* Implement the Feynman shift-add algorithm 
 * The argument z must lie between 1 and 2 inclusive
 */
static void fp_log_shiftadd (mpf_t lg, const mpf_t z, unsigned int prec)
{
	mpf_t zee, ex, tp, su;
	mpf_init (zee);
	mpf_init (tp);
	mpf_init (ex);
	mpf_init (su);
	
	/* Make a copy of the input arguments now! */
	mpf_set (zee, z);
	mpf_set_ui (lg, 0);
	
	mpf_set_ui (tp, 1);
	mpf_set_ui (ex, 1);
	int n;
  	for (n=1; n<3.322*prec; n++)
	{
		// mpf_div_ui (tp, tp, 2);
		mpf_div_2exp (tp, tp, 1);
		while (1)
		{
			mpf_mul (su, ex, tp);
			mpf_add (su, su, ex);
			if (mpf_cmp(su, zee) > 0) break;

			fp_log_2exp (ex, n, prec);
			mpf_add (lg, lg, ex);
			mpf_set (ex, su);
		}
	}
	mpf_clear (su);
	mpf_clear (ex);
	mpf_clear (tp);
	mpf_clear (zee);
}

void fp_log (mpf_t lg, const mpf_t z, unsigned int prec)
{
	mpf_t zee;
	mpf_init (zee);
	mpf_set (zee, z);
	int nexp = 0;

	/* Find power of two by means of bit-shifts */
	while (mpf_cmp_ui(zee, 2) > 0)
	{
		nexp ++;
		mpf_div_ui (zee, zee, 2);
	}

	while (mpf_cmp_ui(zee, 1) < 0)
	{
		nexp --;
		mpf_mul_ui (zee, zee, 2);
	}

#if 1
	/* Apply simple-minded series summation
	 * This appears to be faster than the Feynman algorithm
	 * for the types of numbers and precisions I'm encountering. */
	if (mpf_cmp_d(zee, 1.618) > 0)
	{
		mpf_ui_div (zee, 1, zee);
		mpf_ui_sub (zee, 1, zee);
		fp_log_m1 (lg, zee, prec);
	}
	else
	{
		mpf_ui_sub (zee, 1, zee);
		fp_log_m1 (lg, zee, prec);
		mpf_neg (lg, lg);
	}
#else
	fp_log_shiftadd (lg, zee, prec);
#endif

	/* Add log (2^n) = n log (2) to result. */
	if (0 != nexp)
	{
		fp_log2 (zee, prec);
		if (0 > nexp)
		{
			mpf_mul_ui (zee, zee, -nexp);
			mpf_neg (zee, zee);
		}
		else
		{
			mpf_mul_ui (zee, zee, nexp);
		}
		mpf_add (lg,lg, zee);
	}

	mpf_clear (zee);
}

/**
 * fp_log_ui-- return log(k) for integer k.
 *
 * The values are cached, allowing improved algorithm speeds.
 */
void fp_log_ui (mpf_t lg, unsigned int k, unsigned int prec)
{
	DECLARE_FP_CACHE (log_n);

	if (prec <= fp_one_d_cache_check (&log_n, k))
	{
		fp_one_d_cache_fetch (&log_n, lg, k);
		return;
	}
	
	mpf_set_ui (lg, k);
	fp_log (lg, lg, prec);

	fp_one_d_cache_store (&log_n, lg, k, prec);
}

/* ======================================================================= */
/**
 * cpx_log_m1 -  Floating point logarithm
 * Computes -log(1-z) using Taylor expansion for small z.
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 * Direct cut-n-paste of algos above.
 */

void cpx_log_m1 (cpx_t lg, const cpx_t z, unsigned int prec)
{
	cpx_t zee, z_n, term;

	cpx_init (zee);
	cpx_init (z_n);
	cpx_init (term);

	/* Make copy of argument now! */
	cpx_set (zee, z);
	cpx_mul (z_n, zee, zee);
	cpx_set (lg, zee);
	
	/* Use 10^{-prec} for smallest term in sum */
	mpf_t sqterm, maxterm;
	mpf_init (maxterm);
	mpf_init (sqterm);
	fp_epsilon (maxterm, 2*prec);

	unsigned int n=2;
	while(1)
	{
		cpx_div_ui (term, z_n, n);
		cpx_add (lg, lg, term);
		
		/* don't go no farther than this */
		cpx_mod_sq (sqterm, term);
		if (mpf_cmp (sqterm, maxterm) < 0) break;
		
		n ++;
		cpx_mul (z_n, z_n, zee);
	}
	
	cpx_clear (zee);
	cpx_clear (z_n);
	cpx_clear (term);

	mpf_clear (maxterm);
	mpf_clear (sqterm);
}

void cpx_log (cpx_t lg, const cpx_t z, unsigned int prec)
{
	mpf_t r;
	mpf_init (r);

	/* log (re^{itheta}) = log(r) +  itheta) 
	 * theta = arctan (y/x)
	 */
	cpx_mod_sq (r, z);
	fp_arctan2 (lg[0].im, z[0].im, z[0].re, prec);
	fp_log (lg[0].re, r, prec);

	/* divide by two; using mpf_div_2exp() is faster than mpf_div_ui() */
	// mpf_div_ui (lg[0].re, lg[0].re, 2);
	mpf_div_2exp (lg[0].re, lg[0].re, 1);

	mpf_clear (r);
}

/* ======================================================================= */
/**
 * atan_series -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. 
 * Implements A&S equation 4.4.42 
 * Input value of z should be less than 1/2 for fastest convergence.
 */

static void atan_series (mpf_t atn, const mpf_t zee, unsigned int prec)
{
	mpf_t z_n, zsq, term;

	mpf_init (z_n);
	mpf_init (zsq);
	mpf_init (term);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);

	mpf_mul (zsq, zee, zee);
	mpf_mul (z_n, zee, zsq);
	mpf_set (atn, zee);
		
	unsigned int n=1;
	while(1)
	{
		mpf_div_ui (term, z_n, 2*n+1);
		if (n%2)
		{
			mpf_sub (atn, atn, term);
		}
		else
		{
			mpf_add (atn, atn, term);
		}
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n ++;
		mpf_mul (z_n, z_n, zsq);
	}

	mpf_clear (z_n);
	mpf_clear (zsq);
	mpf_clear (term);

	mpf_clear (maxterm);
}

/**
 * atan2_reduce -- use half-angle reduction to quickly eval arctangent.
 *
 * The half-angle reduction formula is 
 * atan(z) = 2 atan (z/(1+sqrt(1+z^2)))
 *
 * We also get a "free" reduction (no sqrt needed) by using
 * atan (z) = pi/2 - atan (1/z) 
 * whenever z>1.
 *
 * For the half-angle formula to work, its assumed that x>0
 */
static void atan2_reduce (mpf_t atn, const mpf_t y, const mpf_t x, unsigned int prec)
{
	double fy = fabs (mpf_get_d (y));
	double fx = fabs (mpf_get_d (x));
	if (fy > fx)
	{
		mpf_t pih;
		mpf_init (pih);
		fp_pi_half (pih, prec);
		atan2_reduce (atn, x, y, prec);
		mpf_sub (atn, pih, atn);
		mpf_clear (pih);
	}
	if (fy > 0.3*fx)
	{
		mpf_t newx, ysq;
		mpf_init (newx);
		mpf_init (ysq);
		mpf_mul (ysq, y, y);
		mpf_mul (newx, x, x);
		mpf_add (newx, newx, ysq);
		mpf_sqrt (newx, newx);
		mpf_add (newx, newx, x);

		atan2_reduce (atn, y, newx, prec);
		mpf_mul_ui (atn, atn, 2);

		mpf_clear (ysq);
		mpf_clear (newx);
	}
	else
	{
		mpf_t zee;
		mpf_init (zee);

		mpf_div (zee, y, x);
		atan_series (atn, zee, prec);
		mpf_clear (zee);
	}
}

/**
 * fp_arctan -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no serious attempts at optimization. 
 * Implements A&S equation 4.4.42 
 */

void fp_arctan2 (mpf_t atn, const mpf_t y, const mpf_t x, unsigned int prec)
{
	int sgn_x = mpf_sgn(x);
	int sgn_y = mpf_sgn(y);
	
	if (0 < sgn_y)
	{
		if (0 < sgn_x)
		{
			atan2_reduce (atn, y, x, prec);
		}
		else if (0 > sgn_x)
		{
			mpf_t pi, negx;
			mpf_init (pi);
			mpf_init (negx);
			fp_pi (pi, prec);
			mpf_neg (negx, x);

			atan2_reduce (atn, y, negx, prec);
			mpf_sub (atn, pi, atn);

			mpf_clear (pi);
			mpf_clear (negx);
		}
		else
		{
			fp_pi_half (atn, prec);
		}
	}
	else 
	{
		if (0 < sgn_x)
		{
			atan2_reduce (atn, y, x, prec);
		}
		else if (0 > sgn_x)
		{
			mpf_t negpi, negx;
			mpf_init (negpi);
			mpf_init (negx);
			fp_pi (negpi, prec);
			mpf_neg (negpi, negpi);
			mpf_neg (negx, x);

			atan2_reduce (atn, y, negx, prec);
			mpf_sub (atn, negpi, atn);

			mpf_clear (negpi);
			mpf_clear (negx);
		}
		else
		{
			fp_pi_half (atn, prec);
			mpf_neg (atn, atn);
		}
	}
}

void fp_arctan (mpf_t atn, const mpf_t z, unsigned int prec)
{
	mpf_t one;
	mpf_init (one);
	mpf_set_ui (one,1);
	fp_arctan2 (atn, z, one, prec);
	mpf_clear (one);
}

/* ======================================================================= */

void cpx_sqrt (cpx_t rt, const cpx_t z, int prec)
{
	mpf_t modulus;
	mpf_init (modulus);

	cpx_mod_sq (modulus, z);
	mpf_sqrt (modulus, modulus);
	
#ifdef GET_HALF_ANGLE_ARCTAN_STYLE
	mptf_t phase;
	mpf_init (phase);
	fp_arctan2 (phase, z[0].im, z[0].re, prec);
	mpf_div_ui (phase, phase, 2);

	fp_cosine (rt[0].re, phase, prec);
	fp_sine (rt[0].im, phase, prec);
	mpf_clear (phase);
#endif

#define GET_HALF_ANGLE_TRIG_IDENT_STYLE
#ifdef GET_HALF_ANGLE_TRIG_IDENT_STYLE
	/* use half-angle formulas for sine and cosine */
	/* cos A = sqrt(0.5*(1+cos 2A)) */
	cpx_div_mpf (rt, z, modulus);
	mpf_add_ui (rt[0].re, rt[0].re, 1);

	/* divide by two; using mpf_div_2exp() is faster than mpf_div_ui() */
	// mpf_div_ui (rt[0].re, rt[0].re, 2);
	mpf_div_2exp (rt[0].re, rt[0].re, 1);
	mpf_sqrt (rt[0].re, rt[0].re);
	
	/* avoid divide by zero when cosine(half) is zero
	 * i.e. when sine (full angle) is zero */
	if (mpf_cmp_ui(rt[0].im,0))
	{
		/* sin A = sin 2A / (2 cos A) */
		// mpf_div_ui (rt[0].im, rt[0].im, 2);
		mpf_div_2exp (rt[0].im, rt[0].im, 1);
		mpf_div (rt[0].im, rt[0].im, rt[0].re);
	}
	else
	{
		mpf_set_ui (rt[0].im, 1);
	}
#endif

	mpf_sqrt (modulus, modulus);
	cpx_times_mpf (rt, rt, modulus);

	mpf_clear (modulus);
}

/* ======================================================================= */
/*
 * cpx_mpf_pow -- return q^s for complex s, real q.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ... Consider using fp_pow_rc instead, which will 
 * cache values if q and s are heald const, while k is varied.
 */
void cpx_mpf_pow (cpx_t powc, const mpf_t kq, const cpx_t ess, int prec)
{
	mpf_t logkq, mag, pha;
	mpf_init (logkq);
	mpf_init (mag);
	mpf_init (pha);

	/* Domain error is kq is not positive 
	 * We could define this to be exp (i pi s) times result, but I'm lazy */
	if (mpf_cmp_ui (kq, 0) < 0)
	{
		fprintf (stderr, "cpx_mpf_pow() domain error, q<0\n");
		return;
	}

	fp_log (logkq, kq, prec);
	
	/* magnitude is exp(re(s) * log(kq)) */
	mpf_mul (mag, ess->re, logkq);
	
	fp_exp (mag, mag, prec);

	/* phase is im(s) * log(kq)) */
	mpf_mul (pha, ess->im, logkq);

	fp_cosine (powc->re, pha, prec);
	mpf_mul (powc->re, mag, powc->re);
	
	fp_sine (powc->im, pha, prec);
	mpf_mul (powc->im, mag, powc->im);
	
	mpf_clear(logkq);
	mpf_clear(mag);
	mpf_clear(pha);
}

/* ======================================================================= */
/*
 * cpx_pow -- return q^s for complex q, s.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ... If q and s can be held const, while k is varied,
 * then consider using cpx_pow_rc instead, which will cache values 
 * for performance.
 */
void cpx_pow (cpx_t powc, const cpx_t que, const cpx_t ess, int prec)
{
	cpx_t logq;
	cpx_init (logq);

	cpx_log (logq, que, prec);
	cpx_mul (logq, logq, ess);
	cpx_exp (powc, logq, prec);

	cpx_clear(logq);
}

/* ======================================================================= */
/*
 * cpx_pow_ui-- return q^n for complex q, integer n.
 *
 * Uses a log(n) algo, by masking out the bit-string of n.
 * Basically, walk over the bits in n, and multiply by that
 * appropriate power-of-two power for the arg.
 */
void cpx_pow_ui (cpx_t powc, const cpx_t q, unsigned int n)
{
	int k;
	cpx_t qsq;

	cpx_init (qsq);

	cpx_set (qsq, q);
	cpx_set_ui (powc, 1, 0);

	for (k=0; k<32;k++)
	{
		if (0 == n) return;

		if (n & 0x1)
		{
			cpx_mul (powc, powc, qsq);
		}

		n >>= 1;
		cpx_mul(qsq, qsq, qsq);
	}
	
	cpx_clear (qsq);
}

/* ======================================================================= */
/**
 * cpx_ui_pow -- return k^s for complex s, integer k.
 *
 * Uses a brute-force algo: it requires a logarithm, an exp, sin 
 * and cos to be computed, each of which are kinda slow ... 
 */
void cpx_ui_pow (cpx_t powc, unsigned int k, const cpx_t ess, int prec)
{
	mpf_t logkq, mag, pha;
	mpf_init (logkq);
	mpf_init (mag);
	mpf_init (pha);

	fp_log_ui (logkq, k, prec);
	
	/* magnitude is exp(re(s) * log(kq)) */
	mpf_mul (mag, ess->re, logkq);
	
	fp_exp (mag, mag, prec);

	/* phase is im(s) * log(kq)) */
	mpf_mul (pha, ess->im, logkq);

	fp_cosine (powc->re, pha, prec);
	mpf_mul (powc->re, mag, powc->re);
	
	fp_sine (powc->im, pha, prec);
	mpf_mul (powc->im, mag, powc->im);
	
	mpf_clear(logkq);
	mpf_clear(mag);
	mpf_clear(pha);
}

/**
 * cpx_ui_pow_cache -- return k^s for complex s, integer k.
 *
 * If s is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */
void cpx_ui_pow_cache (cpx_t powc, unsigned int k, const cpx_t ess, int prec)
{
	DECLARE_CPX_CACHE (powcache);
	static cpx_t cache_s;
	static int precision = 0;

	if (!precision)
	{
		cpx_init (cache_s);
	}

	if (precision < prec)
	{
		precision = prec;
		cpx_set_prec (cache_s, 3.322*prec+50);
	}

	if (!cpx_eq (ess, cache_s, prec*3.322))
	{
		cpx_one_d_cache_clear (&powcache);
		cpx_set (cache_s, ess);
	}

	if (prec <= cpx_one_d_cache_check (&powcache, k))
	{
		cpx_one_d_cache_fetch (&powcache, powc, k);
		return;
	}
	
	cpx_ui_pow (powc, k, ess, prec);

	cpx_one_d_cache_store (&powcache, powc, k, prec);
}

/* ======================================================================= */
/**
 * fp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q and s is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds. Up to two distinct values
 * of q are cached.
 *
 * Overall, though, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ... 
 */
void fp_pow_rc (cpx_t powc, int k, const mpf_t q, const cpx_t ess, int prec)
{
	DECLARE_CPX_CACHE (powcache_one);
	DECLARE_CPX_CACHE (powcache_two);
	static mpf_t cache_q_one, cache_q_two;
	static cpx_t cache_s_one, cache_s_two;
	static int precision = 0;
	static cpx_cache *next_cache;
	static mpf_t *next_q;
	static cpx_t *next_s;

	if (!precision)
	{
		mpf_init (cache_q_one);
		mpf_init (cache_q_two);
		cpx_init (cache_s_one);
		cpx_init (cache_s_two);
		next_cache = &powcache_one;
		next_q = &cache_q_one;
		next_s = &cache_s_one;
	}

	if (precision < prec)
	{
		precision = prec;
		mpf_set_prec (cache_q_one, 3.322*prec+50);
		mpf_set_prec (cache_q_two, 3.322*prec+50);
		cpx_set_prec (cache_s_one, 3.322*prec+50);
		cpx_set_prec (cache_s_two, 3.322*prec+50);
		cpx_one_d_cache_clear (&powcache_one);
		cpx_one_d_cache_clear (&powcache_two);
	}

	cpx_cache *powcache = NULL;
	if (mpf_eq(q, cache_q_one, prec*3.322) &&
	    cpx_eq(ess, cache_s_one, prec*3.322))
	{
		powcache = &powcache_one;
		if (prec <= cpx_one_d_cache_check (powcache, k))
		{
			cpx_one_d_cache_fetch (powcache, powc, k);
			return;
		}
		next_cache = &powcache_two;
		next_q = &cache_q_two;
		next_s = &cache_s_two;
	}
	
	if (mpf_eq(q, cache_q_two, prec*3.322) &&
	    cpx_eq(ess, cache_s_two, prec*3.322))
	{
		powcache = &powcache_two;
		if (prec <= cpx_one_d_cache_check (powcache, k))
		{
			cpx_one_d_cache_fetch (powcache, powc, k);
			return;
		}
		next_cache = &powcache_one;
		next_q = &cache_q_one;
		next_s = &cache_s_one;
	}
	
	if (NULL == powcache)
	{
		powcache = next_cache;
		cpx_one_d_cache_clear (powcache);
		cpx_one_d_cache_check (powcache, 4);
		mpf_set(*next_q, q);
		cpx_set(*next_s, ess);
	}

	mpf_t kq;
	mpf_init (kq);
	mpf_add_ui (kq, q, k);
	cpx_mpf_pow (powc, kq, ess, prec);
	mpf_clear (kq);

	cpx_one_d_cache_store (powcache, powc, k, prec);
}

void cpx_pow_rc (cpx_t powc, int k, const cpx_t q, const cpx_t ess, int prec)
{
	DECLARE_CPX_CACHE (powcache_one);
	DECLARE_CPX_CACHE (powcache_two);
	static cpx_t cache_q_one, cache_q_two;
	static cpx_t cache_s_one, cache_s_two;
	static int precision = 0;
	static cpx_cache *next_cache;
	static cpx_t *next_q, *next_s;

	if (!precision)
	{
		cpx_init (cache_q_one);
		cpx_init (cache_q_two);
		cpx_init (cache_s_one);
		cpx_init (cache_s_two);
		next_cache = &powcache_one;
		next_q = &cache_q_one;
		next_s = &cache_s_one;
	}

	if (precision < prec)
	{
		precision = prec;
		cpx_set_prec (cache_q_one, 3.322*prec+50);
		cpx_set_prec (cache_q_two, 3.322*prec+50);
		cpx_set_prec (cache_s_one, 3.322*prec+50);
		cpx_set_prec (cache_s_two, 3.322*prec+50);
		cpx_one_d_cache_clear (&powcache_one);
		cpx_one_d_cache_clear (&powcache_two);
	}

	cpx_cache *powcache = NULL;
	if (cpx_eq(q, cache_q_one, prec*3.322) &&
	    cpx_eq(ess, cache_s_one, prec*3.322))
	{
		powcache = &powcache_one;
		if (prec <= cpx_one_d_cache_check (powcache, k))
		{
			cpx_one_d_cache_fetch (powcache, powc, k);
			return;
		}
		next_cache = &powcache_two;
		next_q = &cache_q_two;
		next_s = &cache_s_two;
	}
	
	if (cpx_eq(q, cache_q_two, prec*3.322) &&
	    cpx_eq(ess, cache_s_two, prec*3.322))
	{
		powcache = &powcache_two;
		if (prec <= cpx_one_d_cache_check (powcache, k))
		{
			cpx_one_d_cache_fetch (powcache, powc, k);
			return;
		}
		next_cache = &powcache_one;
		next_q = &cache_q_one;
		next_s = &cache_s_one;
	}
	
	if (NULL == powcache)
	{
		powcache = next_cache;
		cpx_one_d_cache_clear (powcache);
		cpx_one_d_cache_check (powcache, 4);
		cpx_set(*next_q, q);
		cpx_set(*next_s, ess);
	}

	
	cpx_t kq;
	cpx_init (kq);
	mpf_add_ui (kq[0].re, q[0].re, k);
	mpf_set (kq[0].im, q[0].im);
	cpx_pow (powc, kq, ess, prec);
	cpx_clear (kq);

	cpx_one_d_cache_store (powcache, powc, k, prec);
}

/* =============================== END OF FILE =========================== */


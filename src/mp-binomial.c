/*
 * mp-binomial.c
 *
 * High-precison factorials and binomials, using the 
 * Gnu Multiple-precision library.
 * 
 * Copyright (C) 2005, 2006 Linas Vepstas
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
#include "mp-cache.h"
#include "mp-binomial.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* ======================================================================= */
/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */

void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n)
{
	mpz_t acc;
	mpz_init (acc);

	mpz_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpz_mul_ui (acc, poch, i+k);
		mpz_set (poch, acc);
	}

	mpz_clear (acc);
}

/** 
 * i_factorial -- the factorial
 */
// #define USE_LOCAL_FACTORIAL
#ifdef USE_LOCAL_FACTORIAL
void i_factorial (mpz_t fact, unsigned int n)
{
	DECLARE_I_CACHE (cache);

	if (1 >= n)
	{
		mpz_set_ui (fact, 1);
		return;
	}
	int hit = i_one_d_cache_check (&cache, n);
	if (hit)
	{
		i_one_d_cache_fetch (&cache, fact, n);
	}
	else
	{
		i_poch_rising (fact, 1, n);
		i_one_d_cache_store (&cache, fact, n);
	}
}
#else
#define i_factorial mpz_fac_ui
#endif /* USE_LOCAL_FACTORIAL */

/* ====================================================================== */

/* compute and cache inverse factorial */
void fp_inv_factorial (mpf_t inv, unsigned int k, unsigned int prec)
{
	DECLARE_FP_CACHE (infac);

	if (prec <= fp_one_d_cache_check (&infac, k))
	{
		fp_one_d_cache_fetch (&infac, inv, k);
		return;
	}

	mpz_t fac;
	mpz_init (fac);
	mpz_fac_ui (fac, k);
	mpf_set_z (inv, fac);
	mpf_ui_div (inv, 1, inv);
	
	mpz_clear (fac);
	fp_one_d_cache_store (&infac, inv, k, prec);
}

/* ====================================================================== */
/* i_binomial
 * Binomial coefficient (n k)
 */

// #define USE_LOCAL_BINOMIAL
#ifdef USE_LOCAL_BINOMIAL
static void i_binomial_compute (mpz_t bin, unsigned int n, unsigned int k)
{
	mpz_t top, bot;

	if (2*k < n) k = n-k;

	mpz_init (top);
	mpz_init (bot);
	i_poch_rising (top, k+1, n-k);
	i_factorial (bot, n-k); 

	mpz_divexact (bin, top, bot);
	
	mpz_clear (top);
	mpz_clear (bot);
}

static void i_binomial_recurse (mpz_t bin, unsigned int n, unsigned int k)
{
	mpz_t top, bot;

	mpz_init (top);
	mpz_init (bot);

	i_binomial (bot,n-1,k-1);
	i_binomial (top,n-1,k);
	mpz_add (bin, top, bot);
	
	mpz_clear (top);
	mpz_clear (bot);
}

/**
 * i_binomial - return the binomial coefficient
 * Uses a cached value if avalable.
 */ 
void i_binomial (mpz_t bin, unsigned int n, unsigned int k)
{
	DECLARE_I_CACHE (cache);

	if (k > n || 0 > k)
	{
		mpz_set_ui (bin, 0);
		return;
	}

	if (1 >= n)
	{
		mpz_set_ui (bin, 1);
		return;
	}

	if (2*k < n) k = n-k;
	int hit = i_triangle_cache_check (&cache, n, k);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, bin, n, k);
	}
	else
	{
		// i_binomial_compute (bin, n, k);
		i_binomial_recurse (bin, n, k);
		i_triangle_cache_store (&cache, bin, n, k);
	}
}
#else 
#define i_binomial mpz_bin_uiui

#endif /* USE_LOCAL_BINOMIAL */

/* ======================================================================= */

/**
 * i_binomial_sequence -- returns binomial, assumes purely sequeintial access
 * 
 * This routine assumes that the binomial coefficients will be 
 * accessed in an utterly sequential mode, with k running from 
 * zero to n, and n running from zero to k. For sequential access,
 * this routine is very very fast. Otherwise, random access is used
 * which is considerably slower.
 */
void i_binomial_sequence (mpz_t bin, unsigned int n, unsigned int k)
{
	static int curr_n=0, last_k=0;
	DECLARE_I_CACHE (a_cache);
	DECLARE_I_CACHE (b_cache);
	static i_cache *curr_cache = &a_cache;
	static i_cache *next_cache = &b_cache;

	/* Gap in access sequence; fill in the gap */
	if (k > last_k+1 || (n > curr_n && k !=0))
	{
		int j,m;
		if (n == curr_n)
		{
			for (j=last_k+1; j<k; j++)
			{
				i_binomial_sequence (bin, n, j);
			}
		}
		else
		{
			for (j=last_k+1; j<=curr_n; j++)
			{
				i_binomial_sequence (bin, curr_n, j);
			}
			for (m=curr_n+1; m<n; m++)
			{
				for (j=0; j<=m; j++)
				{
					i_binomial_sequence (bin, m, j);
				}
			}
			for (j=0; j<k; j++)
			{
				i_binomial_sequence (bin, n, j);
			}
		}
	}
	
	/* standard sequential access */
	if (k == last_k+1 && n == curr_n)
	{
		/* End of current row -- special case */
		if (k == n)
		{
			mpz_set_ui (bin, 1);
			i_one_d_cache_store (next_cache, bin, k);
			return;
		}

		mpz_t bl;
		mpz_init (bl);

		i_one_d_cache_fetch (curr_cache, bin, k-1);
		i_one_d_cache_fetch (curr_cache, bl, k);
		mpz_add (bin, bin, bl);
		i_one_d_cache_store (next_cache, bin, k);
		mpz_clear (bl);

		last_k = k;
		return;
	}
	
	/* Start a new row */
	if (k == 0 && n == curr_n+1)
	{
		last_k = 0;
		curr_n = n;
		i_cache *tmp = curr_cache;
		curr_cache = next_cache;
		next_cache = tmp;
		i_one_d_cache_check (next_cache, n+1);
		mpz_set_ui (bin, 1);
		i_one_d_cache_store (next_cache, bin, 0);
		return;
	}

	/* Initialize the suential access system. */
	if (0 == n && 0 == k)
	{
		curr_n = 0;
		last_k = 0;
		i_one_d_cache_check (curr_cache, 3);
		i_one_d_cache_check (next_cache, 3);
		mpz_set_ui (bin, 1);
		return;
	}

	/* invalid input */
	if (k > n || 0 > k)
	{
		mpz_set_ui (bin, 0);
		return;
	}

	/* If we got to here, it must be some random access. */
	i_binomial (bin, n, k);
fprintf (stderr, "booooo! n=%d k=%d  currn=%d lastk=%d\n", n,k, curr_n, last_k);
}

/* ======================================================================= */
/* stirling_first - Stirling Numbers of the First kind, 
 * normalized so that they are all positive.
 * Uses dynamically-sized cache.
 */
void i_stirling_first (mpz_t s, unsigned int n, unsigned int k)
{
	DECLARE_I_CACHE (cache);

	/* Trivial case (not in the cache) */
	if (0==k)
	{
		if (0==n) 
		{ 
			mpz_set_ui (s, 1);
		}
		else
		{
			mpz_set_ui (s, 0);
		}
		return;
	}

	if (n<k)
	{
		mpz_set_ui (s, 0);
		return;
	}

	if (n==k)
	{
		mpz_set_ui (s, 1);
		return;
	}

	/* Pull value from cache if it is there */
	int hit = i_triangle_cache_check (&cache, n, k);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, s, n, k);
		return;
	}
	
	/* Use recursion to get new value */
	/* s(n,k) = s(n-1, k-1) + (n-1) * s(n-1, k) */
	unsigned int i;
	mpz_t skm, sk, en;
	mpz_init (skm);
	mpz_init (sk);
	mpz_init (en);
	mpz_set_ui (skm, 0);
	mpz_set_ui (en, n-1);
	for (i=1; i<=n; i++)
	{
		i_stirling_first (sk, n-1, i);
		mpz_mul (s, en, sk);
		mpz_add (s, s, skm);
		i_triangle_cache_store (&cache, s, n, i);
		mpz_set (skm, sk);
	}
	mpz_clear (skm);
	mpz_clear (sk);
	mpz_clear (en);

	i_triangle_cache_fetch (&cache, s, n, k);
}

/* ======================================================================= */
/* A funny off-by-one sum of stirling and binomial */

static void i_stirbin_sum_compute (mpz_t s, unsigned int n, unsigned int m)
{
	unsigned int k;
	
	mpz_t term, stir, bin;
	mpz_init (term);
	mpz_init (stir);
	mpz_init (bin);
	mpz_set_ui (s, 0);
	for (k=m; k<=n; k++)
	{
		i_stirling_first (stir, n, k);
		i_binomial (bin, k,m);
		mpz_mul (term, bin, stir);
		if (k%2)
		{
			mpz_sub (s, s, term);
		}
		else
		{
			mpz_add (s, s, term);
		}
	}
	mpz_clear (term);
	mpz_clear (stir);
	mpz_clear (bin);
}

void i_stirbin_sum (mpz_t s, unsigned int n, unsigned int m)
{
	DECLARE_I_CACHE (cache);

	if (0 >= n)
	{
		mpz_set_ui (s, 1);
		return;
	}

	int hit = i_triangle_cache_check (&cache, n, m);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, s, n, m);
	}
	else
	{
		i_stirbin_sum_compute (s, n, m);
		i_triangle_cache_store (&cache, s, n, m);
	}
}

/* ======================================================================= */
/* stirling_second - Stirling Numbers of the Second kind, 
 * Uses dynamically-sized cache.
 */
void i_stirling_second (mpz_t s, unsigned int n, unsigned int k)
{
	DECLARE_I_CACHE (cache);

	/* Trivial case (not in the cache) */
	if (0==k)
	{
		if (0==n) 
		{ 
			mpz_set_ui (s, 1);
		}
		else
		{
			mpz_set_ui (s, 0);
		}
		return;
	}

	if (n<k)
	{
		mpz_set_ui (s, 0);
		return;
	}

	if (n==k)
	{
		mpz_set_ui (s, 1);
		return;
	}

	/* Pull value from cache if it is there */
	int hit = i_triangle_cache_check (&cache, n, k);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, s, n, k);
		return;
	}
	
	/* Use recursion to get new value */
	/* S(n, k) = S(n-1, k-1) + k * S(n-1, k) */
	unsigned int i;
	mpz_t skm, sk, eye;
	mpz_init (skm);
	mpz_init (sk);
	mpz_init (eye);
	mpz_set_ui (skm, 0);
	for (i=1; i<=n; i++)
	{
		i_stirling_second (sk, n-1, i);
		mpz_set_ui (eye, i);
		mpz_mul (s, eye, sk);
		mpz_add (s, s, skm);
		i_triangle_cache_store (&cache, s, n, i);
		mpz_set (skm, sk);
	}
	mpz_clear (skm);
	mpz_clear (sk);
	mpz_clear (eye);

	i_triangle_cache_fetch (&cache, s, n, k);
}

/* ======================================================================= */
/* binomial transform of power sum */

static void fp_bin_xform_pow_compute (mpf_t bxp, unsigned int n, unsigned int s)
{
	mpz_t bin;
	mpz_init (bin);

	mpf_t vp, term;
	mpf_init (vp);
	mpf_init (term);
	
	mpf_set_ui (bxp, 0);
	unsigned int k;
	for (k=0; k<=n; k++)
	{
		i_binomial (bin, n, k);
		mpf_set_z (term, bin);
		fp_inv_pow (vp, k+1, s);
		mpf_mul (term, term, vp);

		if (k%2)
		{
			mpf_sub (bxp, bxp, term);
		}
		else
		{
			mpf_add (bxp, bxp, term);
		}
	}
	mpz_clear (bin);
	mpf_clear (vp);
	mpf_clear (term);
}

void fp_bin_xform_pow (mpf_t bxp, unsigned int n, unsigned int s)
{
	DECLARE_FP_CACHE (cache);
	if (0 == n)
	{
		mpf_set_ui (bxp, 1); 
		return;
	}
	int hit = fp_triangle_cache_check (&cache, n+s, s);
	if (hit)
	{
		fp_triangle_cache_fetch (&cache, bxp, n+s, s);
	}
	else
	{
		fp_bin_xform_pow_compute (bxp, n, s);
		fp_triangle_cache_store (&cache, bxp, n+s, s, 1);
	}
}

/* ======================================================================= */
/** 
 * fp_harmonic -- The harmonic number
 */
void fp_harmonic (mpf_t harm, unsigned int n, unsigned int prec)
{
	DECLARE_FP_CACHE (cache);

	if (1 >= n)
	{
		mpf_set_ui (harm, 1);
		return;
	}
	int have_prec = fp_one_d_cache_check (&cache, n);
	if (prec <= have_prec)
	{
		fp_one_d_cache_fetch (&cache, harm, n);
		return;
	}
	
	unsigned int istart = n-1;
	have_prec = fp_one_d_cache_check (&cache, istart);
	while ((have_prec < prec) && (1 < istart))
	{
		istart--;
		have_prec = fp_one_d_cache_check (&cache, istart);
	}

	unsigned int i;
	fp_harmonic (harm, istart, prec);

	mpf_t term;
	mpf_init (term);
	for (i=istart+1; i<=n; i++)
	{
		mpf_set_ui (term, 1);
		mpf_div_ui (term, term, i);
		mpf_add (harm, harm, term);
		fp_one_d_cache_store (&cache, harm, i, prec);
	}
	mpf_clear (term);
}

/* ======================================================================= */
/* fp_poch_rising
 * rising pochhammer symbol (x)_n, for real values of x and integer n.
 *
 * Brute force, simple.
 */

void fp_poch_rising_d (mpf_t poch, double x, unsigned int n)
{
	mpf_t term;
	mpf_init (term);

	mpf_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpf_set_d (term, x+i);
		mpf_mul (poch, poch, term);
	}

	mpf_clear (term);
}

void fp_poch_rising (mpf_t poch, mpf_t x, unsigned int n)
{
	mpf_t term;
	mpf_init (term);

	/* Make copy of arg NOW! */
	mpf_set (term,x);

	mpf_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpf_mul (poch, poch, term);
		mpf_add_ui (term, term, 1);
	}

	mpf_clear (term);
}

/* ======================================================================= */
/* cpx_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */

void cpx_poch_rising_d (cpx_t poch, double re_s, double im_s, unsigned int n)
{
	cpx_t term, acc;
	cpx_init (term);
	cpx_init (acc);

	cpx_set_ui (acc, 1, 0);
	
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpf_set_d (term[0].re, re_s+i);
		mpf_set_d (term[0].im, im_s);

		cpx_mul (acc, acc, term);
	}

	cpx_set (poch, acc);
	cpx_clear (term);
	cpx_clear (acc);
}

/* ======================================================================= */
/* cpx_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */

void cpx_poch_rising (cpx_t poch, const cpx_t ess, unsigned int n)
{
	/* Unlikely special case */
	if (0 == n)
	{
		cpx_set_ui (poch, 1, 0);
		return;
	}
	
	cpx_t term, acc;
	cpx_init (term);
	cpx_init (acc);

	cpx_set (acc, ess);
	cpx_set (term, ess);

	unsigned int i;
	for (i=1; i<n; i++)
	{
		mpf_add_ui (term[0].re, term[0].re, 1);
		cpx_mul (acc, acc, term);
	}

	cpx_set (poch, acc);
	cpx_clear (term);
	cpx_clear (acc);
}

/* ======================================================================= */
/* fp_binomial
 * Binomial coefficient 
 */

void fp_binomial_d (mpf_t bin, double s, unsigned int k)
{
	mpf_t top, bot;
	mpz_t fac;

	mpf_init (top);
	mpf_init (bot);
	mpz_init (fac);
	fp_poch_rising_d (top, s-k+1, k);
	i_factorial (fac, k); 
	mpf_set_z (bot, fac);

	mpf_div (bin, top, bot);
	
	mpf_clear (top);
	mpf_clear (bot);
	mpz_clear (fac);
}

/* ======================================================================== */
/* cpx_binomial
 * Complex binomial coefficient
 */

void cpx_binomial_d (cpx_t bin, double re_s, double im_s, unsigned int k)
{
	cpx_t top;
	cpx_init (top);
	
	cpx_poch_rising_d (top, re_s-k+1, im_s, k);
	
	mpz_t ifac;
	mpz_init (ifac);
	i_factorial (ifac, k); 
	
	mpf_t fac;
	mpf_init (fac);
	mpf_set_z (fac, ifac);
	mpz_clear (ifac);

	cpx_div_mpf (bin, top, fac);
	
	mpf_clear (fac);
	cpx_clear (top);
}

void cpx_binomial (cpx_t bin, const cpx_t ess, unsigned int k)
{
	if (0 >= k) {
		cpx_set_ui (bin, 1, 0);
		return;
	}
			  
	cpx_t top, bot;
	cpx_init (top);
	cpx_init (bot);
	
	cpx_set (bot, ess);
	mpf_sub_ui (bot[0].re, bot[0].re, k-1);

	cpx_poch_rising (top, bot, k);
	
	mpz_t ifac;
	mpz_init (ifac);
	i_factorial (ifac, k); 

	mpf_t fac;
	mpf_init(fac);
	mpf_set_z (fac, ifac);

	cpx_div_mpf (bin, top, fac);
	
	mpf_clear (fac);
	mpz_clear (ifac);
	
	cpx_clear (top);
	cpx_clear (bot);
}

void cpx_binomial_sum_cache (cpx_t bin, const cpx_t ess, unsigned int k)
{
	DECLARE_CPX_CACHE (cache);
	static int cache_bits = 0;
	static cpx_t cache_s;
	
	if (!cache_bits)
	{
		cpx_init (cache_s);
		cpx_set_ui (cache_s, 1, 0);
	}
	int prec = mpf_get_default_prec();
	if (cache_bits < prec)
	{
		cpx_set_prec (cache_s, prec);
		cache_bits = prec;
	}

	/* First, check if this is the same s value as before */
	if (!cpx_eq (cache_s, ess, prec))
	{
		cpx_one_d_cache_clear (&cache);
		cpx_set (cache_s, ess);
	}
	
	/* Check the local cache */
	int have_prec = cpx_one_d_cache_check (&cache, k);
	if (have_prec >= prec)
	{
		cpx_one_d_cache_fetch (&cache, bin, k);
		return;
	}

	cpx_t sn;
	cpx_init (sn);
	cpx_add_ui (sn, ess, k, 0);
	cpx_binomial (bin, sn, k);
	cpx_one_d_cache_store (&cache, bin, k, prec);
	cpx_clear (sn);
}


/* =============================== END OF FILE =========================== */


/*
 * mp-consts.c
 *
 * High-precison constants, using the 
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
#include "mp-binomial.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* ======================================================================= */
/**
 * fp_half_sqrt_three - return 0.5*sqrt(3)= 0.86602...
 */

void fp_half_sqrt_three (mpf_t sqt)
{
	static unsigned int init=0;
	static mpf_t cached_sqt;

	if (init)
	{
		mpf_set (sqt, cached_sqt);
		return;
	}
	mpf_init (cached_sqt);

	mpf_set_ui (sqt, 3);
	mpf_sqrt (sqt, sqt);
	mpf_div_ui (sqt, sqt, 2);
	mpf_set (cached_sqt, sqt);

	init =1;
}

/* ======================================================================= */
/**
 * fp_e - return e=2.718281828...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */

extern void fp_exp_helper (mpf_t ex, const mpf_t z, unsigned int prec);

void fp_e (mpf_t e, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_e;

	if (precision >= prec)
	{
		mpf_set (e, cached_e);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_e);
	}

	mpf_set_prec (cached_e, 3.322*prec +50);

	mpf_t one;
	mpf_init (one);
	mpf_set_ui (one, 1);
	fp_exp_helper (cached_e, one, prec);
	mpf_set (e, cached_e);

	mpf_clear (one);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_pi - return pi=3.14159... 
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force Machin formula
 */
void fp_pi (mpf_t pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_pi;

	if (precision >= prec)
	{
		mpf_set (pi, cached_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_pi);
	}
	mpf_set_prec (cached_pi, 3.322*prec +50);

	/* Simple-minded Machin formula */
	mpf_t tmp;
	mpf_init(tmp);
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 5);
	fp_arctan (pi, tmp, prec);

	mpf_mul_ui (pi, pi, 4);

	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 239);
	fp_arctan (tmp, tmp, prec);

	mpf_sub (pi, pi, tmp);
	mpf_mul_ui (pi, pi, 4);
	mpf_clear (tmp);

	mpf_set (cached_pi, pi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_two_pi - return 2pi = 2 * 3.14159... 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_two_pi (mpf_t two_pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_two_pi;

	if (precision >= prec)
	{
		mpf_set (two_pi, cached_two_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_two_pi);
	}
	mpf_set_prec (cached_two_pi, 3.322*prec +50);

	fp_pi (two_pi, prec);
	mpf_mul_ui (two_pi, two_pi, 2);
	mpf_set (cached_two_pi, two_pi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_two_over_pi - return 2/pi = 2 / 3.14159... 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_two_over_pi (mpf_t two_over_pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_two_over_pi;

	if (precision >= prec)
	{
		mpf_set (two_over_pi, cached_two_over_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_two_over_pi);
	}
	mpf_set_prec (cached_two_over_pi, 3.322*prec +50);

	fp_pi (two_over_pi, prec);
	mpf_ui_div (two_over_pi, 2, two_over_pi);
	mpf_set (cached_two_over_pi, two_over_pi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_pi_half - return pi/2 = 0.5 * 3.14159... 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_pi_half (mpf_t pih, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_pih;

	if (precision >= prec)
	{
		mpf_set (pih, cached_pih);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_pih);
	}
	mpf_set_prec (cached_pih, 3.322*prec +50);

	fp_pi (pih, prec);
	mpf_div_ui (pih, pih, 2);
	mpf_set (cached_pih, pih);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_sqrt_two_pi - return sqrt(2pi) = sqrt (2 * 3.14159...) 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_sqrt_two_pi (mpf_t sqtpi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_sqtpi;

	if (precision >= prec)
	{
		mpf_set (sqtpi, cached_sqtpi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_sqtpi);
	}
	mpf_set_prec (cached_sqtpi, 3.322*prec +50);

	fp_two_pi (sqtpi, prec);
	mpf_sqrt (sqtpi, sqtpi);
	mpf_set (cached_sqtpi, sqtpi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_log_two_pi - return log(2pi) = log(2 * 3.14159...) 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_log_two_pi (mpf_t ltp, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_ltp;

	if (precision >= prec)
	{
		mpf_set (ltp, cached_ltp);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_ltp);
	}
	mpf_set_prec (cached_ltp, 3.322*prec +50);

	fp_two_pi (ltp, prec);
	fp_log (ltp, ltp, prec);
	mpf_set (cached_ltp, ltp);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_log2 - return log(2)=0.69...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */

void fp_log2 (mpf_t l2, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_log2;

	if (precision >= prec)
	{
		mpf_set (l2, cached_log2);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_log2);
	}
	mpf_set_prec (cached_log2, 3.322*prec +50);

	mpf_t two;
	mpf_init (two);
	mpf_set_ui (two, 2);
	fp_log (cached_log2, two, prec);
	mpf_set (l2, cached_log2);

	mpf_clear (two);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_e_pi - return e^pi 
 * @prec - number of decimal places of precision
 *
 * Uses simple, low-brow formula
 */
void fp_e_pi (mpf_t e_pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_e_pi;

	if (precision >= prec)
	{
		mpf_set (e_pi, cached_e_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_e_pi);
	}
	mpf_set_prec (cached_e_pi, 3.322*prec +50);

	fp_pi (e_pi, prec);
	fp_exp (e_pi, e_pi, prec);

	mpf_set (cached_e_pi, e_pi);
	precision = prec;
}


/* ======================================================================= */
/** 
 * fp_euler - return Euler-Mascheroni const
 * @prec - number of decimal places of precision
 *
 */
static void fp_euler_mascheroni_limit (mpf_t gam, unsigned int n, unsigned int prec)
{
	mpf_t maxterm;
	mpf_init (maxterm);
	mpf_set_ui (maxterm, 1);

	mpf_t z_n, twon, term, tmp, fact;
	mpf_init (z_n);
	mpf_init (twon);
	mpf_init (term);
	mpf_init (tmp);
	mpf_init (fact);
	mpf_set_ui (twon, 1);
	mpf_mul_2exp(twon, twon, n);
	mpf_mul (z_n, twon, twon);
	mpf_set_ui (fact, 1);
	mpf_div_ui (fact, fact, 2);
	mpf_set (gam, twon);

	/* The k=1 term is handled above in init */
	unsigned int k=2;
	while(1)
	{
		fp_harmonic (tmp, k, prec);
		mpf_mul (term, z_n, tmp);
		mpf_mul (term, term, fact);

		mpf_add (gam, gam, term);

		/* XXX in fact, we can terminate this sum a lot earlier,
		 * if we wanted to, by comparing the sum to the number of
		 * desired digits ... */
		/* Don't go no farther than this */
		if (mpf_cmp (term, maxterm) < 0) break;

		k ++;
		mpf_mul (z_n, z_n, twon);
		mpf_div_ui (fact, fact, k);
	}

	fp_exp (tmp, twon, prec);
	mpf_div (gam, gam, tmp);

	fp_log2 (tmp, prec);
	mpf_mul_ui (tmp, tmp, n);
	mpf_sub (gam, gam, tmp);
	
	mpf_clear (z_n);
	mpf_clear (twon);
	mpf_clear (term);
	mpf_clear (tmp);
	mpf_clear (fact);

	mpf_clear (maxterm);
}

static void fp_euler_mascheroni_compute (mpf_t gam, unsigned int prec)
{
	/* power value, goes as log log n */
	// double en = log (prec*log(10)) / log (2.0);
	double en = log (prec*3.322);
	// en = 1.442695041 * (en - log (en));
	en = 1.442695041 * en;
	int n = (int) (en+1.0);
	fp_euler_mascheroni_limit (gam, n, prec);
}

void fp_euler_mascheroni (mpf_t gam, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_gam;

	if (precision >= prec)
	{
		mpf_set (gam, cached_gam);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_gam);
	}
	mpf_set_prec (cached_gam, 3.322*prec +50);

	fp_euler_mascheroni_compute (gam, prec);
	mpf_set (cached_gam, gam);
	precision = prec;
}

/* ======================================================================= */
/** 
 * fp_zeta_half - return zeta (1/2)
 * @prec - number of decimal places of precision
 *
 */
static void fp_zeta_half_compute (mpf_t gam, unsigned int prec)
{
	cpx_t ess, zeta;
	cpx_init (ess);
	cpx_init (zeta);
	
	mpf_set_d (ess[0].re, 0.5);
	mpf_set_ui (ess[0].im, 0);

	cpx_borwein_zeta (zeta, ess, prec);
	mpf_set (gam, zeta[0].re);

	cpx_clear (ess);
	cpx_clear (zeta);
}

void fp_zeta_half (mpf_t gam, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_gam;

	if (precision >= prec)
	{
		mpf_set (gam, cached_gam);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_gam);
	}
	mpf_set_prec (cached_gam, 3.322*prec +50);

	fp_zeta_half_compute (gam, prec);
	mpf_set (cached_gam, gam);
	precision = prec;
}

/* =============================== END OF FILE =========================== */


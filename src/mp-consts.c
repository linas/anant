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

#include <pthread.h>
#include <math.h>
#include <stdbool.h>
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
// multi-threading locks.
// All the constants share one lock, there should be no contention.

static pthread_spinlock_t mp_const_lock;
static pthread_spinlock_t mp_pi_lock;
static pthread_spinlock_t mp_euler_lock;
static pthread_spinlock_t mp_zeta_lock;
__attribute__((constructor)) void fp_e_ctor(void)
{
	pthread_spin_init(&mp_const_lock, PTHREAD_PROCESS_PRIVATE);
	pthread_spin_init(&mp_pi_lock, PTHREAD_PROCESS_PRIVATE);
	pthread_spin_init(&mp_euler_lock, PTHREAD_PROCESS_PRIVATE);
	pthread_spin_init(&mp_zeta_lock, PTHREAD_PROCESS_PRIVATE);
}

/* ======================================================================= */
/**
 * fp_half_sqrt_three - return 0.5*sqrt(3)= 0.86602...
 */
void fp_half_sqrt_three (mpf_t sqt, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_sqt;

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (sqt, cached_sqt);
		pthread_spin_unlock(&mp_const_lock);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_sqt);
	}

	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;
	mpf_set_prec (cached_sqt, bits);

	mpf_set_ui (sqt, 3);
	mpf_sqrt (sqt, sqt);
	mpf_div_ui (sqt, sqt, 2);
	mpf_set (cached_sqt, sqt);

	precision = prec;
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (e, cached_e);
		pthread_spin_unlock(&mp_const_lock);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_e);
	}

	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;
	mpf_set_prec (cached_e, bits);

	mpf_t one;
	mpf_init2 (one, bits);
	mpf_set_ui (one, 1);
	fp_exp_helper (cached_e, one, prec);
	mpf_set (e, cached_e);

	mpf_clear (one);
	precision = prec;
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_pi_lock);
	if (precision >= prec)
	{
		mpf_set (pi, cached_pi);
		pthread_spin_unlock(&mp_pi_lock);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_pi);
	}
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;
	mpf_set_prec (cached_pi, bits);

	/* Simple-minded Machin formula */
	mpf_t tmp;
	mpf_init2 (tmp, bits);
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
	pthread_spin_unlock(&mp_pi_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (two_pi, cached_two_pi);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (two_over_pi, cached_two_over_pi);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (pih, cached_pih);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (sqtpi, cached_sqtpi);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (ltp, cached_ltp);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (l2, cached_log2);
		pthread_spin_unlock(&mp_const_lock);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_log2);
	}
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;
	mpf_set_prec (cached_log2, bits);

	mpf_t two;
	mpf_init2 (two, bits);
	mpf_set_ui (two, 2);
	fp_log (cached_log2, two, prec);
	mpf_set (l2, cached_log2);

	mpf_clear (two);
	precision = prec;
	pthread_spin_unlock(&mp_const_lock);
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

	pthread_spin_lock(&mp_const_lock);
	if (precision >= prec)
	{
		mpf_set (e_pi, cached_e_pi);
		pthread_spin_unlock(&mp_const_lock);
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
	pthread_spin_unlock(&mp_const_lock);
}


/* ======================================================================= */
/**
 * fp_euler - return Euler-Mascheroni const
 * @prec - number of decimal places of precision
 *
 */
static void fp_euler_mascheroni_limit (mpf_t gam, unsigned int n, unsigned int prec)
{
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;
	mpf_t maxterm;
	mpf_init2 (maxterm, bits);
	mpf_set_ui (maxterm, 1);

	mpf_t z_n, twon, term, tmp, fact;
	mpf_init2 (z_n, bits);
	mpf_init2 (twon, bits);
	mpf_init2 (term, bits);
	mpf_init2 (tmp, bits);
	mpf_init2 (fact, bits);
	mpf_set_ui (twon, 1);
	mpf_mul_2exp(twon, twon, n);
	mpf_mul (z_n, twon, twon);
	mpf_set_ui (fact, 1);
	mpf_div_ui (fact, fact, 2);
	mpf_set (gam, twon);

	/* The k=1 term is handled above in init */
	unsigned int k=2;
	while (1)
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

	pthread_spin_lock(&mp_euler_lock);
	if (precision >= prec)
	{
		mpf_set (gam, cached_gam);
		pthread_spin_unlock(&mp_euler_lock);
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
	pthread_spin_unlock(&mp_euler_lock);
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

	pthread_spin_lock(&mp_zeta_lock);
	if (precision >= prec)
	{
		mpf_set (gam, cached_gam);
		pthread_spin_unlock(&mp_zeta_lock);
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
	pthread_spin_unlock(&mp_zeta_lock);
}

/* =============================== END OF FILE =========================== */

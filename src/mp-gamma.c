
/*
 * mp-gamma.c
 *
 * Compute gamma function for various complex arguments
 *
 * Copyright (C) 2006 Linas Vepstas
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

#include <gmp.h>

#include "mp-binomial.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* ================================================= */
/*
 * fp_lngamma -- compute log of gamma for real argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 * Slightly modified: 
 * ln Gamma (z+2)=z(1-gamma)+sum_{n=2}^\infty (zeta(n)-1) z^n/n 
 * Valid input must have 0 < x < 4
 * and ideally 1.5 < x < 2.5 for fastest convergence 
 */
static void reduced_lngamma (mpf_t gam, const mpf_t ex, int prec)
{
	int n;
	mpf_t z, zn, term;

	mpf_init (z);
	mpf_init (zn);
	mpf_init (term);

	/* make copy of input argument now! */
	mpf_set (z, ex);

	mpf_set_ui (gam, 0);
	
	mpf_sub_ui (z, z, 2);
	mpf_mul (zn, z,z);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);
	
	n=2;
	while (1)
	{
		fp_zeta (term, n, prec);
		mpf_sub_ui (term, term, 1);
		mpf_mul (term, term, zn);
		mpf_div_ui (term, term, n);
		if (n%2)
		{
			mpf_sub (gam, gam, term);
		}
		else
		{
			mpf_add (gam, gam, term);
		}

		/* don't go no farther than this */
		mpf_abs(term, term);
		if (mpf_cmp (term, maxterm) < 0) break;

		mpf_mul (zn,zn,z);
		n++;
	}

	fp_euler_mascheroni (term, prec);
	mpf_sub_ui (term, term, 1);
	mpf_mul (term, term, z);
	mpf_sub (gam, gam, term);

	mpf_clear (z);
	mpf_clear (zn);
	mpf_clear (term);
}

/* ================================================= */
/*
 * cpx_lngamma -- compute log of gamma for complex argument
 *
 * Same code as above, extended for complex arguments
 * Uses simple, quickly converging algo-- A&S 6.1.33
 * Slightly modified: 
 * ln Gamma (z+2)=z(1-gamma)+sum_{n=2}^\infty (zeta(n)-1) z^n/n 
 * Valid input must have |ex-2| < 2
 * and ideally |ex-2| < 0.5 for fastest convergence 
 */
static void cpx_reduced_lngamma (cpx_t gam, const cpx_t ex, int prec)
{
	int n;
	cpx_t z, zn, term;

	cpx_init (z);
	cpx_init (zn);
	cpx_init (term);

	/* make copy of input argument now! */
	cpx_set (z, ex);

	cpx_set_ui (gam, 0, 0);
	
	cpx_sub_ui (z, z, 2, 0);
	cpx_mul (zn, z, z);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, 2*prec);
	
	n=2;
	while (1)
	{
		mpf_set_ui (term[0].im, 0);
		fp_zeta (term[0].re, n, prec);
		mpf_sub_ui (term[0].re, term[0].re, 1);
		mpf_div_ui (term[0].re, term[0].re, n);
		cpx_mul (term, term, zn);
		
		if (n%2)
		{
			cpx_sub (gam, gam, term);
		}
		else
		{
			cpx_add (gam, gam, term);
		}

		/* don't go no farther than this */
		cpx_mod_sq(term[0].re, term);
		if (mpf_cmp (term[0].re, maxterm) < 0) break;

		cpx_mul (zn,zn,z);
		n++;
	}

	mpf_set_ui (term[0].im, 0);
	fp_euler_mascheroni (term[0].re, prec);
	mpf_sub_ui (term[0].re, term[0].re, 1);
	cpx_mul (term, term, z);
	cpx_sub (gam, gam, term);

	cpx_clear (z);
	cpx_clear (zn);
	cpx_clear (term);
}

/* ================================================= */
/* 
 * gamma function, valid only for 0 < x < 4 
 * Otherwise explodes...
 */ 
static inline void reduced_gamma (mpf_t gam, const mpf_t ex, int prec)
{
	reduced_lngamma (gam, ex, prec);
	fp_exp (gam, gam, prec);
}

/* ================================================= */
/* 
 * fp_gamma
 * Use pochhammer to get into range of 1.5 < z < 2.5,
 * Then use the reduced summation formula. The range
 * of 1.5 < z < 2.5 is particularly quickly convergent
 */
void fp_gamma (mpf_t gam, const mpf_t z, int prec)
{
	mpf_t zee;
	mpf_init (zee);

	/* make a copy of the input arg NOW! */
	mpf_set (zee, z);
	
	/* double-presision used, this code doesn't need to 
	 * be all that accurate; just need a reasonable int part. 
	 */
	double flo = mpf_get_d (zee);
	if (flo > 2.5)
	{
		unsigned int intpart = (unsigned int) floor (flo-1.0);
		
		/* The goal of the next if statement is to make sure that
		 * -0.5<zee<0.5, which helps maintain excellent convergence.
		 */
		if (flo-intpart < 1.5) intpart --;
		mpf_sub_ui (zee, zee, intpart);
		fp_poch_rising (gam, zee, intpart);
	}
	else if (flo < 1.5)
	{
		unsigned int intpart = (unsigned int) floor (2.0-flo);

		/* The goal of the next if statement is to make sure that
		 * -0.5<zee<0.5, which helps maintain excellent convergence.
		 */
		if (flo+intpart < 1.5) intpart ++;
		fp_poch_rising (gam, zee, intpart);
		mpf_ui_div (gam, 1, gam);

		mpf_add_ui (zee, zee, intpart);
	}
	else
	{
		mpf_set_ui (gam, 1);
	}

	mpf_t rgamma;
	mpf_init (rgamma);
	reduced_gamma (rgamma, zee, prec);
	
	mpf_mul (gam, gam, rgamma);

	mpf_clear (zee);
	mpf_clear (rgamma);
}

void mpf_gamma_cache (mpf_t gam, const mpf_t z, int prec)
{
	static mpf_t cache_z, cache_gam;
	static int precision = 0;
	int redo = 0;

	if (!precision)
	{
		mpf_init (cache_z);
		mpf_init (cache_gam);
	}
	if (precision < prec)
	{
		mpf_set_prec (cache_z, 3.322*prec+50);
		mpf_set_prec (cache_gam, 3.322*prec+50);
		precision = prec;
		redo = 1;
	}

	if (redo || !mpf_eq (z, cache_z, prec*3.322))
	{
		mpf_set (cache_z, z);
		fp_gamma (gam, z, prec);
		mpf_set (cache_gam, gam);
	}
	else
	{
		mpf_set (gam, cache_gam);
	}
}

/* ================================================= */
/* 
 * gamma function, but valid only for -2 < ImZ < 2 
 * Use pochhammer to get into range of 1.5 < Re Z < 2.5,
 * Then use the reduced summation formula. The range
 * of 1.5 < ReZ < 2.5 is particularly quickly convergent
 */ 

static void cpx_reduced_gamma (cpx_t gam, const cpx_t z, int prec)
{
	cpx_t zee;
	cpx_init (zee);

	/* make a copy of the input arg NOW! */
	cpx_set (zee, z);
	
	/* double-presision used, this code doesn't need to 
	 * be all that accurate; just need a reasonable int part. 
	 */
	double flo = mpf_get_d (zee[0].re);
	if (flo > 2.5)
	{
		unsigned int intpart = (unsigned int) floor (flo-1.0);
		
		/* The goal of the next if statement is to make sure that
		 * -0.5<zee<0.5, which helps maintain excellent convergence.
		 */
		if (flo-intpart < 1.5) intpart --;
		cpx_sub_ui (zee, zee, intpart, 0);
		cpx_poch_rising (gam, zee, intpart);
	}
	else if (flo < 1.5)
	{
		unsigned int intpart = (unsigned int) floor (2.0-flo);

		/* The goal of the next if statement is to make sure that
		 * -0.5<zee<0.5, which helps maintain excellent convergence.
		 */
		if (flo+intpart < 1.5) intpart ++;
		cpx_poch_rising (gam, zee, intpart);
		cpx_recip (gam, gam);

		cpx_add_ui (zee, zee, intpart, 0);
	}
	else
	{
		cpx_set_ui (gam, 1, 0);
	}

	cpx_t rgamma;
	cpx_init (rgamma);

	cpx_reduced_lngamma (rgamma, zee, prec);
	cpx_exp (rgamma, rgamma, prec);
	cpx_mul (gam, gam, rgamma);

	cpx_clear (zee);
	cpx_clear (rgamma);
}

/* ================================================= */
/* 
 * gamma function for general complex argument
 * Use the multiplication theorem
 */ 
void cpx_gamma (cpx_t gam, const cpx_t z, int prec)
{
	/* step one: find out how big the imaginary part is */
	double img = fabs(mpf_get_d (z[0].im));
	int m = (int) (img + 1.0);

	if (1 == m)
	{
		cpx_reduced_gamma (gam, z, prec);
		return;
	}
	
	/* step two: prepare to use the multiplicatin theorem */
	cpx_t zee, mzee, term, acc;
	cpx_init (zee);
	cpx_init (mzee);
	cpx_init (term);
	cpx_init (acc);

	/* Copy the input arg NOW! */
	cpx_set (mzee, z);
	
	mpf_t frac;
	mpf_init (frac);

	cpx_div_ui (zee, mzee, m);

	/* frac = 1/m */
	mpf_set_ui (frac, 1);
	mpf_div_ui (frac, frac, m);
	
	cpx_set_ui (acc, 1, 0);

	/* Apply the multiplication theorem */
	int k;
	for (k=0; k<m; k++)
	{
		cpx_reduced_gamma (term, zee, prec);
		cpx_mul (acc, acc, term);
		mpf_add (zee[0].re, zee[0].re, frac);
	}

	/* Multiply by scaling factors */
	fp_pi (frac, prec);
	mpf_mul_ui (frac, frac, 2);
	mpf_sqrt (frac, frac);   /* XXX could avoid sqrt if m is odd ... */
	mpf_pow_ui (frac, frac, m-1);
	cpx_div_mpf (acc, acc, frac);

	/* mz - 0.5 */
	mpf_set_ui (frac, 1);
	// mpf_div_ui (frac, frac, 2);
	mpf_div_2exp (frac, frac, 1);  /* frac = frac/2, which is faster than mpf_div_ui(frac,frac,2); */
	mpf_sub (mzee[0].re, mzee[0].re,frac);

	/* m^(mz-0.5) */
	mpf_set_ui (frac, m);
	cpx_mpf_pow (term, frac, mzee, prec);
	cpx_mul (acc, acc, term);

	cpx_set (gam, acc);

	cpx_clear (zee); 
	cpx_clear (mzee); 
	cpx_clear (term); 
	cpx_clear (acc); 
	mpf_clear (frac);
}

void cpx_gamma_cache (cpx_t gam, const cpx_t z, int prec)
{
	static cpx_t cache_z, cache_gam;
	static int precision = 0;
	int redo = 0;

	if (!precision)
	{
		cpx_init (cache_z);
		cpx_init (cache_gam);
	}
	
	if (precision < prec)
	{
		cpx_set_prec (cache_z, 3.322*prec+50);
		cpx_set_prec (cache_gam, 3.322*prec+50);
		precision = prec;
		redo = 1;
	}

	if (redo || !cpx_eq (z, cache_z, prec*3.322))
	{
		cpx_set (cache_z, z);
		cpx_gamma (gam, z, prec);
		cpx_set (cache_gam, gam);
	}
	else
	{
		cpx_set (gam, cache_gam);
	}
}

/* ==================  END OF FILE ===================== */

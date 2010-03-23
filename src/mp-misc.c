/*
 * mp-misc.c
 *
 * High-precison misc functions, using the 
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2005,2006,2007 Linas Vepstas
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
#include "mp-complex.h"
#include "mp-misc.h"

void i_prt (const char * str, mpz_t val)
{
	printf ("%s", str);
	mpz_out_str (stdout, 10, val);
}

void fp_prt (const char * str, mpf_t val)
{
	printf ("%s", str);
	mpf_out_str (stdout, 10, 60, val);
}

void cpx_prt (const char * str, const cpx_t val)
{
	printf ("%s", str);
	mpf_out_str (stdout, 10, 30, val[0].re);
	printf (" + i ");
	mpf_out_str (stdout, 10, 30, val[0].im);
}

void ecpx_prt (const char * str, const cpx_t val)
{
	fprintf (stderr, "%s", str);
	mpf_out_str (stderr, 10, 30, val[0].re);
	fprintf (stderr, " + i ");
	mpf_out_str (stderr, 10, 30, val[0].im);
}

/* ===================================================== */
/**
 * fp_epsilon - return 10^{-prec} 
 */
void fp_epsilon (mpf_t eps, int prec)
{
	static int cache_prec = -1;
	static mpf_t cache_eps;

	if (-1 == cache_prec)
	{
		mpf_init (cache_eps);
	}

	if (prec == cache_prec)
	{
		mpf_set (eps, cache_eps);
		return;
	}
	if (cache_prec < prec)
	{
		mpf_set_prec (cache_eps, 3.322*prec+50);
	}

	/* double mex = ((double) prec) * log (10.0) / log(2.0); */
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t one;
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (cache_eps, one, imax);

	mpf_set (eps, cache_eps);
	cache_prec = prec;
	mpf_clear (one);
}

/* ===================================================== */

/* prec is the decimal precison (number of decimal places) */
/* nterms is the number of an's to compute */
void set_bits (int prec, int nterms)
{
	/* Compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* The variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */
	int bits = (int) (v + 300 + 3*nterms);

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
}

/* ===================================================== */

int get_prec (cpx_t epsi, unsigned int prec)
{
	long rex, imx;
	mpf_get_d_2exp (&rex, epsi[0].re);
	mpf_get_d_2exp (&imx, epsi[0].im);
	rex = -0.30103 *rex;
	imx = -0.30103 *imx;
	if (imx && imx < rex) rex = imx;
	if (0 == rex) rex = imx;
	if (0 == rex) rex = prec; 
	if (mpf_cmp_d (epsi[0].re, 0.1) > 0) rex = 0;
	if (mpf_cmp_d (epsi[0].re, -0.1) < 0) rex = 0;
	
	if (mpf_cmp_d (epsi[0].im, 0.1) > 0) rex = 0;
	if (mpf_cmp_d (epsi[0].im, -0.1) < 0) rex = 0;
	
	return rex;
}

/* Print number of digits by which value differes from 
 * previous call to this routine.
 */

int last_change(const cpx_t curr, unsigned int prec)
{
	static cpx_t prev;
	static int init = 0;

	if (!init)
	{
		init = 1;
		cpx_init (prev);
	}

	/* Set the precision (number of binary bits) */
	int nbits = 3.322*prec+5;
	cpx_set_prec (prev, nbits);

	cpx_sub (prev, prev, curr);

	printf ("prec=%d ", prec);

	long rex = get_prec (prev, prec);
	printf ("change=%ld\n", rex);

	cpx_set (prev, curr);
	return rex;
}

/* =============================== END OF FILE =========================== */


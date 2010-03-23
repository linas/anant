/*
 * mp-hyper.c
 *
 * High-precison Hypergeometric functions, using the 
 * Gnu Multiple-precision library.
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
#include <stdlib.h>

#include <gmp.h>
#include "mp-complex.h"
#include "mp-misc.h"

/* ======================================================================= */
/**
 * cpx_confluent -- Confluent hypergeometric function
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void 
cpx_confluent (cpx_t em, cpx_t a, cpx_t b, cpx_t z, unsigned int prec)
{
	mpf_t fact;
	cpx_t zee, z_n, term;
	cpx_t ay, be, poch_a, poch_b;

	mpf_init (fact);

	cpx_init (zee);
	cpx_init (z_n);
	cpx_init (term);
	cpx_init (ay);
	cpx_init (be);
	cpx_init (poch_a);
	cpx_init (poch_b);

	/* Make copy of arguments now! */
	cpx_set (zee, z);
	cpx_set (z_n, zee);
	
	cpx_set (ay, a);
	cpx_set (poch_a, a);
	cpx_set (be, b);
	cpx_set (poch_b, b);
	
	cpx_set_ui (em, 1, 0);
	mpf_set_ui (fact, 1);

	mpf_t mag;
	mpf_init (mag);
	
	/* Use 10^{-2 * prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);
	mpf_mul (maxterm, maxterm, maxterm);

	unsigned int n=1;
	while(1)
	{
		cpx_div_mpf (term, z_n, fact);
		cpx_mul (term, term, poch_a);
		cpx_div (term, term, poch_b);
		cpx_add (em, em, term);
		
		/* Don't go no farther than this */
		cpx_mod_sq (mag, term);
		if (mpf_cmp (mag, maxterm) < 0) break;
		
		n++;
		cpx_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);

		mpf_add_ui (ay[0].re, ay[0].re, 1);
		mpf_add_ui (be[0].re, be[0].re, 1);
		cpx_mul (poch_a, poch_a, ay);
		cpx_mul (poch_b, poch_b, be);
	}
	
	mpf_clear (fact);
	mpf_clear (mag);
	mpf_clear (maxterm);
	
	cpx_clear (zee);
	cpx_clear (z_n);
	cpx_clear (term);
	cpx_clear (ay);
	cpx_clear (be);
	cpx_clear (poch_a);
	cpx_clear (poch_b);
}

/* =============================== END OF FILE =========================== */


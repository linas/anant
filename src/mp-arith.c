/*
 * mp-arith.c
 *
 * High-precison number-theoretic arithmetic series, using the
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2016 Linas Vepstas
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
#include "mp-misc.h"

/* ======================================================================= */
/*
 * sigma_one
 * Sum of divisors. That is, sigma_one(n) = sum_{d|n} d
 * See https://en.wikipedia.org/wiki/Divisor_function
 * Uses cached values.
 *
 * Brute force, very simple. There is a simple and much faster recursive
 * algorithm for this, if we have a list of prime numbers on hand.
 * However, the below is adequate for current purposes (mostly because
 * practical work with the partition function below never gets past
 * a value of about n=50K or 100K or so, so the brute-force sum is tenable).
 */
static void sigma_one_z_nocache (mpz_t sum, unsigned int n)
{
	mpz_set_ui (sum, n);
	unsigned int d;
	unsigned int ns = n / 2;
	for (d=1; d<=ns; d++)
	{
		if (n%d) continue;
		mpz_add_ui (sum, sum, d);
	}
}

void sigma_one_z (mpz_t sum, unsigned int n)
{
	DECLARE_I_CACHE(sigone);

	if (i_one_d_cache_check(&sigone, n))
	{
		mpz_init (sum);
		i_one_d_cache_fetch(&sigone, sum, n);
		return;
	}

	sigma_one_z_nocache(sum, n);
	i_one_d_cache_store(&sigone, sum, n);
}

// ===========================================================
/*
 * partition function.
 * See https://en.wikipedia.org/wiki/Partition_(number_theory)
 * Uses cached values.
 *
 * The problem here is that the partition function overflows
 * a 64-bit int around n=400, and a 128-bit int around n=1400.
 *
 * Brute force, simple.
 */
void partition_z (mpz_t sum, unsigned int n)
{
	DECLARE_I_CACHE(parti);

	if (0 == n)
	{
		mpz_set_ui(sum, 1);
		return;
	}

	if (i_one_d_cache_check(&parti, n))
	{
		i_one_d_cache_fetch(&parti, sum, n);
		return;
	}

	mpz_t sig; mpz_init (sig);
	mpz_t part; mpz_init(part);
	mpz_t term; mpz_init(term);

	mpz_set_ui (sum, 0);
	unsigned int k;
	for (k=0; k<n; k++)
	{
		sigma_one_z(sig, n-k);
		partition_z(part, k);
		mpz_mul(term, sig, part);
		mpz_add (sum, sum, term);
	}
	mpz_div_ui(sum, sum, n);
	mpz_clear(sig);
	mpz_clear(part);
	mpz_clear(term);

	i_one_d_cache_store(&parti, sum, n);
}

// ===========================================================
// #define RUN_TEST
#ifdef RUN_TEST

int main()
{
	mpz_t part; mpz_init (part);

	for (int n=1; n< 20; n++)
	{
		partition_z(part, n);
		unsigned long p = mpz_get_ui(part);
		printf("n=%d p=%ld\n", n, p);
	}
}

#endif

/* =============================== END OF FILE =========================== */

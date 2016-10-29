/*
 * mp-arith.h
 *
 * High-precison partition function, and other number-theoretic
 * arithmetic series, using the Gnu Multiple-precision library.
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

#include <gmp.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * sigma_one
 * Sum of divisors. That is, sigma_one(n) = sum_{d|n} d
 * See https://en.wikipedia.org/wiki/Divisor_function
 * Uses cached values.
 *
 * Brute force, simple.
 */
void sigma_one_z (mpz_t poch, unsigned int n);

/**
 * partition function.
 * See https://en.wikipedia.org/wiki/Partition_(number_theory)
 * Uses cached values.
 *
 * The problem here is that the partition function overflows
 * a 64-bit int around n=400, and a 128-bit int around n=1600.
 *
 * Brute force, simple.
 */
void partition_z (mpz_t poch, unsigned int n);

#ifdef  __cplusplus
};
#endif


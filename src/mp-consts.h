/*
 * mp-consts.h
 *
 * High-precison constants, using the 
 * Gnu Multiple-precision library.
 * Uses caching so as to speed up subsequent derefs.
 *
 * Copyright (C) 2005 Linas Vepstas
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
 * fp_half_sqrt_three - return sqrt(3)/2 = 0.86602...
 */
void fp_half_sqrt_three (mpf_t sqt);

/**
 * fp_e - return e=2.718281828...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */
void fp_e (mpf_t e, unsigned int prec);

/**
 * fp_pi - return pi=3.14159... 
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force Machin formula
 */
void fp_pi (mpf_t pi, unsigned int prec);
void fp_two_pi (mpf_t pi, unsigned int prec);
void fp_pi_half (mpf_t pihalf, unsigned int prec);
void fp_sqrt_two_pi (mpf_t sqtpi, unsigned int prec);
void fp_log_two_pi (mpf_t ltpi, unsigned int prec);
void fp_two_over_pi (mpf_t tpi, unsigned int prec);

/**
 * fp_e_pi - return e^pi 
 * @prec - number of decimal places of precision
 *
 * Uses simple, low-brow formula
 */
void fp_e_pi (mpf_t e_pi, unsigned int prec);

/**
 * fp_log2 - return log(2)=0.693147181...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */
void fp_log2 (mpf_t l2, unsigned int prec);

/**
 * fp_euler - return Euler-Mascheroni const
 * @prec - number of decimal places of precision
 *
 */
void fp_euler_mascheroni (mpf_t gam, unsigned int prec);

/**
 * fp_zeta_half - return zeta(1/2)
 * @prec - number of decimal places of precision
 *
 * Current algorithm is terribly slow.
 */
void fp_zeta_half (mpf_t zeta, unsigned int prec);

#ifdef  __cplusplus
};
#endif


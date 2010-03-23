/*
 * mp-gamma.h
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

#include <gmp.h>
#include "mp-complex.h"

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * fp_gamma -- compute Gamma(x)=factorial(x-1) for real argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 *
 * The caching version skips the calculation, if called again with
 * the same value of ex (up to the nprec precision bits).
 */
void fp_gamma (mpf_t gam, const mpf_t ex, int prec);
void fp_gamma_cache (mpf_t gam, const mpf_t ex, int prec);

/**
 * cpx_gamma -- compute Gamma(x)=factorial(x-1) for complex argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 *
 * The caching version skips the calculation, if called again with
 * the same value of ex (up to the nprec precision bits).
 */
void cpx_gamma (cpx_t gam, const cpx_t ex, int prec);
void cpx_gamma_cache (cpx_t gam, const cpx_t ex, int prec);

#ifdef  __cplusplus
};
#endif

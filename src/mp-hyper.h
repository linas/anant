/*
 * mp-hyper.h
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

#include <gmp.h>
#include "mp-complex.h"

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * cpx_confluent -- Confluent hypergeometric function
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void cpx_confluent (cpx_t em, cpx_t a, cpx_t b, cpx_t z, unsigned int prec);

#ifdef  __cplusplus
};
#endif


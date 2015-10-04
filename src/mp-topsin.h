/*
 * mp-topsin.h
 *
 * High-precison Topologists Sine function coefficients
 * using the Gnu Multiple-precision library.
 *
 * Copyright (C) 2014 Linas Vepstas
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
 *
 */

#include <gmp.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * Series coeffcients of the topologists sine function.  This computes
 * and returns the coefficients a_k from the series
 *
 *    sin(2pi / (1+x)) = sum_k=0^\infty a_k x^k
 *
 * The singularity is placed at x=-1 because that is it's natural
 * location for number-theoretic applications.
 */
void topsin_series (mpf_t a_k, unsigned int k, unsigned int prec);


#ifdef  __cplusplus
};
#endif

/* =============================== END OF FILE =========================== */

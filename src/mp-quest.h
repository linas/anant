/*
 * mp-quest.h
 *
 * High-precison Minkowski Question mark, Stern-Brocot Tree, etc.
 * using the Gnu Multiple-precision library.
 *
 * Copyright (C) 2010 Linas Vepstas
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
 * Minkowski Question Mark function. Valid argument range
 * is 0 =< x =< 1. 
 */
void question_mark (mpf_t qmark, const mpf_t x, unsigned int prec);

/**
 * Inverse of the Minkowski Question Mark function. Valid argument range
 * is 0 =< x =< 1. Implemented with a simple, fast bit-counting
 * algorithm.
 */
void question_inverse (mpf_t qinv, const mpf_t x, unsigned int prec);

#ifdef  __cplusplus
};
#endif

/* =============================== END OF FILE =========================== */

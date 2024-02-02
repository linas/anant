/**
 * mp-zeroiso.h
 *
 * Isolate complex zeros of a polynomial.
 *
 * Copyright (C) 2024 Linas Vepstas
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

#ifndef __MP_ZEROISO_H__
#define __MP_ZEROISO_H__

#ifdef  __cplusplus
extern "C" {
#endif

#include "mp-complex.h"

/* =============================================== */
/**
 * cpx_isolate_roots.
 * Isolate all zeros of a complex-valued square-free polynomial.
 * Returns a list of disks, such that each disk contains exactly
 * one zero (one root).
 *
 * Implements
 *    Michael Sagraloff, Chee K. Yap, "A Simple But Exact and Efficient
 *    Algorithm for Complex Root Isolation" (2011)
 *    https://cs.nyu.edu/exact/doc/complex.pdf
 *
 * Requirements:
 * 1) The polynomial must be square-free, i.e. it cannot have
 *    degenerate zeros, i.e. two zero located at the same place.
 * 2) A means of evaluating the derivative, to all orders, must be
 *    provided.
 *
 * @poly polynomial whose zeros are to be found.
 *       poly takes z and k as input, returns f^(k)(z) as output.
 * @degree Degree of the polynomial.
 * @boxll, @boxur bounding box lower-left and upper-right coordinates.
 *       The search for zeros will be performed inside this box.
 *
 * Fills in:
 * @centers Array of centers of disks, each disk containing one zero.
 * @radii   Array of radii of the disks.
 * Both of these must be provided by the caller, and must be of length
 * at least equal to degree of the polynomial.
 *
 * Returns number of zeros found.
 */
int cpx_isolate_roots(
              void (*poly)(cpx_t f, int deriv, cpx_t z, void* args),
              int degree,
              cpx_t boxll, cpx_t boxur,
              cpx_t* centers, mpf_t* radii,
              void* args);

#ifdef  __cplusplus
};
#endif

#endif /* __MP_ZEROISO_H__ */

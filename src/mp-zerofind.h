/**
 * mp-zerofind.h
 *
 * Locate zeros of a function. 
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
 */

#ifndef __MP_ZEROFIND_H__
#define __MP_ZEROFIND_H__

#ifdef  __cplusplus
extern "C" {
#endif

#include "mp-complex.h"

/* =============================================== */
/**
 * cpx_find_zero.
 * Numerically locate the zero of a complex-valued function.
 * 
 * @func function whose zeros are to be found.
 *       func takes z as input, returns f as output. 
 *       'nprec' is the suggested decimal precisiton at which 'fun'
 *       should perform its calculations.
 * @initial_z initial suggestion for the location of the zero.
 * @e1, @e2 initial suggestions for a bounding ellipse. These are
 *       taken to be two vectors, specifying the major and minor axes of
 *       an ellipse, centered at 'initial_z'. The true zero is presumed
 *       to lie inside of, or at least, close to, this initial ellipse.
 * @ndigits number of decimal digits of accuracy to which the zero
 *       should be searched for.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 *
 * @returns 0 if result is valid, else an error code.
 *
 * This implements Powell's method, slightly adapted; the adaptations
 * are meant to imporve convergence when 'func' is extremely noisy, 
 * i.e. when any sort of quardatic behaviour is obscured by
 * high-freqency noise.
 *
 * A description of Powell's method can be found in Press, Teukolsky,
 * Vetterling, Flannery, "Numerical Recipes in C, 2nd ed.", Cambridge 
 * U Press, 1999.
 */
int cpx_find_zero(cpx_t result,
              void (*func)(cpx_t f, cpx_t z, int nprec),
              cpx_t initial_z,
              cpx_t e1, cpx_t e2,
              int ndigits, int nprec);

#ifdef  __cplusplus
};
#endif

#endif /* __MP_ZEROFIND_H__ */

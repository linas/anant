/*
 * mp-misc.h
 *
 * High-precison misc functions, using the 
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2005,2006,2007 Linas Vepstas
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

void i_prt (const char * str, mpz_t val);
void fp_prt (const char * str, mpf_t val);
void cpx_prt (const char * str, const cpx_t val);
void ecpx_prt (const char * str, const cpx_t val);

/**
 * fp_epsilon - return 10^{-prec} 
 *
 * Intended for use as a max-precision iteration stopper
 * Value is cached.
 */
void fp_epsilon (mpf_t eps, int prec);

/* prec is the decimal precison (number of decimal places) */
/* nterms is the number of an's to compute */
void set_bits (int prec, int nterms);

/**
 * Print number of digits by which value differs from 
 * previous call to this routine.
 */

int last_change(const cpx_t curr, unsigned int prec);

/**
 * get_prec -- return log_10 (epsi) approximately rounded
 */
int get_prec (cpx_t epsi, unsigned int prec);

#ifdef  __cplusplus
};
#endif


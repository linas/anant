/* 
 * mp-complex.h
 *
 * Utility routines for handling complex values
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

#ifndef __MP_COMPLEX_H__
#define __MP_COMPLEX_H__

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct {
	mpf_t re;
	mpf_t im;
} __cpx_struct;

typedef __cpx_struct cpx_t[1];

static inline void cpx_init (cpx_t z)
{
	mpf_init (z[0].re);
	mpf_init (z[0].im);
}

static inline void cpx_init2 (cpx_t z, mp_bitcnt_t bits)
{
	mpf_init2 (z[0].re, bits);
	mpf_init2 (z[0].im, bits);
}

static inline void cpx_clear (cpx_t z)
{
	mpf_clear (z[0].re);
	mpf_clear (z[0].im);
}

static inline void cpx_set (cpx_t z, const cpx_t y)
{
	mpf_set (z[0].re, y[0].re);
	mpf_set (z[0].im, y[0].im);
}

static inline void cpx_set_ui (cpx_t z, unsigned long x, unsigned long y)
{
	mpf_set_ui (z[0].re, x);
	mpf_set_ui (z[0].im, y);
}

static inline void cpx_set_d (cpx_t z, double x, double y)
{
	mpf_set_d (z[0].re, x);
	mpf_set_d (z[0].im, y);
}

static inline void cpx_set_mpf (cpx_t z, const mpf_t x, const mpf_t y)
{
	mpf_set (z[0].re, x);
	mpf_set (z[0].im, y);
}

static inline void cpx_add (cpx_t sum, const cpx_t a, const cpx_t b)
{
	mpf_add (sum[0].re, a[0].re, b[0].re);
	mpf_add (sum[0].im, a[0].im, b[0].im);
}

static inline void cpx_add_d (cpx_t sum, const cpx_t a, double rb, double ib)
{
	mpf_t tmp;
	mpf_init (tmp);
	mpf_set_d (tmp, rb);
	mpf_add (sum[0].re, a[0].re, tmp);
	mpf_set_d (tmp, ib);
	mpf_add (sum[0].im, a[0].im, tmp);
	mpf_clear (tmp);
}

static inline void cpx_add_ui (cpx_t sum, const cpx_t a, unsigned long rb, unsigned long ib)
{
	if (rb) mpf_add_ui (sum[0].re, a[0].re, rb);
	if (ib) mpf_add_ui (sum[0].im, a[0].im, ib);
}

static inline void cpx_add_mpf (cpx_t sum, const cpx_t a, const mpf_t b)
{
	mpf_add (sum[0].re, a[0].re, b);
	mpf_set (sum[0].im, a[0].im);
}

static inline void cpx_sub (cpx_t dif, const cpx_t a, const cpx_t b)
{
	mpf_sub (dif[0].re, a[0].re, b[0].re);
	mpf_sub (dif[0].im, a[0].im, b[0].im);
}

static inline void cpx_sub_ui (cpx_t sum, const cpx_t a, unsigned long rb, unsigned long ib)
{
	if (rb) mpf_sub_ui (sum[0].re, a[0].re, rb);
	if (ib) mpf_sub_ui (sum[0].im, a[0].im, ib);
}

static inline void cpx_ui_sub (cpx_t sum, unsigned long ra, unsigned long ia, const cpx_t b)
{
	if (ra) mpf_ui_sub (sum[0].re, ra, b[0].re);
	if (ia) mpf_ui_sub (sum[0].im, ia, b[0].im);
}

static inline void cpx_sub_mpf (cpx_t sum, const cpx_t a, const mpf_t b)
{
	mpf_sub (sum[0].re, a[0].re, b);
	mpf_set (sum[0].im, a[0].im);
}

static inline void cpx_neg (cpx_t neg, const cpx_t a)
{
	mpf_neg (neg[0].re, a[0].re);
	mpf_neg (neg[0].im, a[0].im);
}

static inline void cpx_conj (cpx_t neg, const cpx_t a)
{
	mpf_neg (neg[0].im, a[0].im);
}

/**
 * cpx_mul -- prod = a * b
 */
static inline void cpx_mul (cpx_t prod, const cpx_t a, const cpx_t b)
{
	mp_bitcnt_t bits = mpf_get_prec(prod[0].re) + 8;
	mpf_t pre, pim, tmp;
	mpf_init2 (pre, bits);
	mpf_init2 (pim, bits);
	mpf_init2 (tmp, bits);
	
	mpf_mul (tmp, a[0].im, b[0].im);
	mpf_mul (pre, a[0].re, b[0].re);
	mpf_sub (pre, pre, tmp);
	
	mpf_mul (tmp, a[0].im, b[0].re);
	mpf_mul (pim, a[0].re, b[0].im);
	mpf_add (pim, pim, tmp);
	
	mpf_set (prod[0].re, pre);
	mpf_set (prod[0].im, pim);

	mpf_clear (pre);
	mpf_clear (pim);
	mpf_clear (tmp);
}

/**
 * cpx_times_i -- z = a*i
 */
static inline void cpx_times_i (cpx_t z, const cpx_t a)
{
	mp_bitcnt_t bits = mpf_get_prec(z[0].re) + 8;
	mpf_t tmp;
	mpf_init2 (tmp, bits);
	mpf_set (tmp, a[0].re);
	mpf_neg (z[0].re, a[0].im);
	mpf_set (z[0].im, tmp);
	mpf_clear (tmp);
}

/**
 * cpx_times_mpf -- prod = a * b
 */
static inline void cpx_times_mpf (cpx_t prod, const cpx_t a, const mpf_t b)
{
	mpf_mul (prod[0].re, a[0].re, b);
	mpf_mul (prod[0].im, a[0].im, b);
}

/**
 * cpx_times_ui -- prod = a * b
 */
static inline void cpx_times_ui (cpx_t prod, const cpx_t a, unsigned long b)
{
	mpf_mul_ui (prod[0].re, a[0].re, b);
	mpf_mul_ui (prod[0].im, a[0].im, b);
}

/**
 * cpx_times_d -- prod = a * b
 */
static inline void cpx_times_d (cpx_t prod, const cpx_t a, double b)
{
	mpf_t tmp;
	mpf_init (tmp);
	mpf_set_d (tmp, b);
	cpx_times_mpf(prod, a, tmp);
	mpf_clear (tmp);
}

/**
 * cpx_recip -- recip = 1/z
 */
static inline void cpx_recip (cpx_t recip, const cpx_t z)
{
	mp_bitcnt_t bits = mpf_get_prec(recip[0].re) + 8;
	mpf_t mag,tmp;
	mpf_init2 (mag, bits);
	mpf_init2 (tmp, bits);
	mpf_mul (mag, z[0].re, z[0].re);
	mpf_mul (tmp, z[0].im, z[0].im);
	mpf_add (mag, mag, tmp);
	mpf_ui_div (mag, 1, mag);
	mpf_mul (recip[0].re, z[0].re, mag);
	mpf_mul (recip[0].im, z[0].im, mag);
	mpf_neg (recip[0].im, recip[0].im);
	
	mpf_clear (mag);
	mpf_clear (tmp);
}

/**
 * cpx_div -- ratio = a/b
 */
static inline void cpx_div (cpx_t ratio, const cpx_t a, const cpx_t b)
{
	mp_bitcnt_t bits = mpf_get_prec(ratio[0].re) + 8;
	cpx_t recip;
	cpx_init2 (recip, bits);

	cpx_recip (recip, b);
	cpx_mul (ratio, a, recip);
	cpx_clear (recip);
}

/**
 * cpx_div_mpf -- ratio = a/b
 */
static inline void cpx_div_mpf (cpx_t ratio, const cpx_t a, const mpf_t b)
{
	mpf_div (ratio[0].re, a[0].re, b);
	mpf_div (ratio[0].im, a[0].im, b);
}

static inline void cpx_div_ui (cpx_t ratio, const cpx_t a, unsigned long b)
{
	mpf_div_ui (ratio[0].re, a[0].re, b);
	mpf_div_ui (ratio[0].im, a[0].im, b);
}

/**
 * Divide a by 2^n place result in ratio
 */
static inline void cpx_div_2exp (cpx_t ratio, const cpx_t a, unsigned long n)
{
	mpf_div_2exp (ratio[0].re, a[0].re, n);
	mpf_div_2exp (ratio[0].im, a[0].im, n);
}

/**
 * cpx_mod_sq -- modulus squared of z
 */
static inline void cpx_mod_sq (mpf_t mod, const cpx_t z)
{
	mp_bitcnt_t bits = mpf_get_prec(mod) + 8;
	mpf_t tmp;
	mpf_init2 (tmp, bits);
	
	mpf_mul (tmp, z[0].im, z[0].im);
	mpf_mul (mod, z[0].re, z[0].re);
	mpf_add (mod, mod, tmp);
	mpf_clear (tmp);
}

/**
 * cpx_abs -- Absolute value of z (modulus of z)
 */
static inline void cpx_abs (mpf_t mod, const cpx_t z)
{
	cpx_mod_sq (mod, z);
	mpf_sqrt (mod, mod);
}

/**
 * Return true if the first nbits of both the imaginary and real parts are equal
 */
static inline int cpx_eq (const cpx_t a, const cpx_t b, unsigned int nbits)
{
	return (mpf_eq(a[0].re, b[0].re, nbits)) && (mpf_eq(a[0].im, b[0].im, nbits));
}

/**
 * Get the real, imaginary parts, as double-precision floats
 */
static inline double cpx_get_re (const cpx_t z)
{
	return mpf_get_d(z[0].re);
}

static inline double cpx_get_im (const cpx_t z)
{
	return mpf_get_d(z[0].im);
}


static inline void cpx_set_prec (cpx_t a, unsigned int nbits)
{
	mpf_set_prec (a[0].re, nbits);
	mpf_set_prec (a[0].im, nbits);
}

#ifdef  __cplusplus
};
#endif

#endif /* __MP_COMPLEX_H__ */

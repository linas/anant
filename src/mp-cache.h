/*
 * mp-cache.h
 *
 * Array cachine functions for arrays of high precision quantities
 * Gnu Multiple-precision library.
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
#include "mp-complex.h"

/* ======================================================================= */
/* Cache management */

typedef struct {
	unsigned int nmax;
	mpz_t *cache;
	char *ticky;
	short disabled;
} i_cache;


#define DECLARE_I_CACHE(name)         \
	static i_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = 0}

/** i_one_d_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int i_one_d_cache_check (i_cache *c, unsigned int n);

/** 
 * i_one_d_cache_fetch - fetch value from cache
 */
static inline void i_one_d_cache_fetch (i_cache *c, mpz_t val, unsigned int n)
{
	if (c->disabled) return;
	mpz_set (val, c->cache[n]);
}

/**
 * i_one_d_cache_store - store value in cache
 */
static inline void i_one_d_cache_store (i_cache *c, const mpz_t val, unsigned int n)
{
	if (c->disabled) return;
	mpz_set (c->cache[n], val);
	c->ticky[n] = 1;
}

void i_one_d_cache_clear (i_cache *c);

/* ======================================================================= */
/** i_triangle_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int i_triangle_cache_check (i_cache *c, unsigned int n, unsigned int k);

/** 
 * i_triangle_cache_fetch - fetch value from cache
 */
static inline void i_triangle_cache_fetch (i_cache *c, mpz_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpz_set (val, c->cache[idx+k]);
}


/**
 * i_triangle_cache_store - store value in cache
 */
static inline void i_triangle_cache_store (i_cache *c, const mpz_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpz_set (c->cache[idx+k], val);
	c->ticky[idx+k] = 1;
}

/* ======================================================================= */
/* Cache management */
/* pure cut-n-paste of he integer variant */

typedef struct {
	unsigned int nmax;
	mpq_t *cache;
	char *ticky;
} q_cache;

#define DECLARE_Q_CACHE(name)         \
	static q_cache name = {.nmax=0, .cache=NULL, .ticky=NULL}

/** q_one_d_cache_check() -- check if mpq_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int q_one_d_cache_check (q_cache *c, unsigned int n);

/** 
 * q_one_d_cache_fetch - fetch value from cache
 */
static inline void q_one_d_cache_fetch (q_cache *c, mpq_t val, unsigned int n)
{
	mpq_set (val, c->cache[n]);
}

/**
 * q_one_d_cache_store - store value in cache
 */
static inline void q_one_d_cache_store (q_cache *c, const mpq_t val, unsigned int n)
{
	mpq_set (c->cache[n], val);
	c->ticky[n] = 1;
}

/* ======================================================================= */
/* Cache management */
/* Almost a cut-n-paste of above, but using fp instead */

typedef struct {
	unsigned int nmax;
	mpf_t *cache;
	int *precision; /* base-10 precision of cached value */
} fp_cache;


#define DECLARE_FP_CACHE(name)         \
	static fp_cache name = {.nmax=0, .cache=NULL, .precision=NULL}

/** fp_one_d_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int fp_one_d_cache_check (fp_cache *c, unsigned int n);

/** 
 * fp_one_d_cache_fetch - fetch value from cache
 */
static inline void fp_one_d_cache_fetch (fp_cache *c, mpf_t val, unsigned int n)
{
	mpf_set (val, c->cache[n]);
}

/**
 * fp_one_d_cache_store - store value in cache
 */
static inline void fp_one_d_cache_store (fp_cache *c, const mpf_t val, unsigned int n, int prec)
{
	mpf_set_prec (c->cache[n], 3.22*prec+50);
	mpf_set (c->cache[n], val);
	c->precision[n] = prec;
}

void fp_one_d_cache_clear (fp_cache *c);

/* ======================================================================= */
/** fp_triangle_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int fp_triangle_cache_check (fp_cache *c, unsigned int n, unsigned int k);

/** 
 * fp_triangle_cache_fetch - fetch value from cache
 */
static inline void fp_triangle_cache_fetch (fp_cache *c, mpf_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpf_set (val, c->cache[idx+k]);
}

/**
 * fp_triangle_cache_store - store value in cache
 */
static inline void fp_triangle_cache_store (fp_cache *c, const mpf_t val, 
					 unsigned int n, unsigned int k, int prec)
{
	unsigned int idx = n * (n+1) /2 ;
	mpf_set_prec (c->cache[idx+k], 3.22*prec+50);
	mpf_set (c->cache[idx+k], val);
	c->precision[idx+k] = prec;
}

/* ======================================================================= */
/* Cache management */
/* A cut-n-paste of above, but using cpx instead */

typedef struct {
	unsigned int nmax;
	cpx_t *cache;
	int *precision; /* base-10 precision of cached value */
} cpx_cache;


#define DECLARE_CPX_CACHE(name)         \
	static cpx_cache name = {.nmax=0, .cache=NULL, .precision=NULL}

/** cpx_one_d_cache_check() -- check if cpx_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int cpx_one_d_cache_check (cpx_cache *c, unsigned int n);

/** 
 * cpx_one_d_cache_fetch - fetch value from cache
 */
static inline void cpx_one_d_cache_fetch (cpx_cache *c, cpx_t val, unsigned int n)
{
	cpx_set (val, c->cache[n]);
}

/**
 * cpx_one_d_cache_store - store value in cache
 */
static inline void cpx_one_d_cache_store (cpx_cache *c, const cpx_t val, unsigned int n, int prec)
{
	cpx_set_prec (c->cache[n], 3.22*prec+50);
	cpx_set (c->cache[n], val);
	c->precision[n] = prec;
}

void cpx_one_d_cache_clear (cpx_cache *c);

/* =============================== END OF FILE =========================== */


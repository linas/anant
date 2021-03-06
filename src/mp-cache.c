/*
 * mp-cache.c
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include "mp-cache.h"

/* ======================================================================= */
/* Cache management */

/** i_one_d_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int i_one_d_cache_check (i_cache *c, unsigned int n)
{
	if (c->disabled) return 0;
	pthread_spin_lock(&c->lock);

	if ((n <= c->nmax) && 0 != c->nmax)
	{
		int gotit = c->ticky[n];
		pthread_spin_unlock(&c->lock);
		return gotit;
	}

	unsigned int newsize = 1.5*n+2;
	mpz_t* new_cache = (mpz_t *) malloc (newsize * sizeof (mpz_t));
	if (c->nmax) memcpy(new_cache, c->cache, (c->nmax+1) * sizeof(mpz_t));

	char* new_ticky = (char *) malloc (newsize * sizeof(char));
	if (c->nmax) memcpy(new_ticky, c->ticky, (c->nmax+1) * sizeof(char));

	unsigned int en;
	unsigned int nstart = c->nmax+1;
	if (0 == c->nmax) nstart = 0;
	for (en = nstart; en < newsize; en++)
	{
		mpz_init (new_cache[en]);
		new_ticky[en] = 0;
	}

	/* Now swap out the old and new */
	mpz_t* old_cache = c->cache;
	c->cache = new_cache;
	char* old_ticky = c->ticky;
	c->ticky = new_ticky;
	c->nmax = newsize-1;
	pthread_spin_unlock(&c->lock);

	if (old_cache) free(old_cache);
	if (old_ticky) free(old_ticky);
	return 0;
}

/* ======================================================================= */
/** i_triangle_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int i_triangle_cache_check (i_cache *c, unsigned int n, unsigned int k)
{
	if (c->disabled) return 0;
	pthread_spin_lock(&c->lock);

	if ((n <= c->nmax) && 0 != c->nmax)
	{
		unsigned int idx = n * (n+1) /2 ;
		int gotit = c->ticky[idx+k];
		pthread_spin_unlock(&c->lock);
		return gotit;
	}

	if (0 == n) n = 1;
	unsigned int newsize = (n+1)*(n+2)/2;
	mpz_t* new_cache = (mpz_t *) malloc (newsize * sizeof (mpz_t));
	if (c->nmax) memcpy(new_cache, c->cache, (c->nmax+1) * sizeof(mpz_t));

	char* new_ticky = (char *) malloc (newsize * sizeof(char));
	if (c->nmax) memcpy(new_ticky, c->ticky, (c->nmax+1) * sizeof(char));

	unsigned int en;
	unsigned int nstart = c->nmax + 1;
	if (0 == c->nmax) nstart = 0;
	for (en = nstart; en <= n; en++)
	{
		unsigned int j;
		unsigned int idx = en * (en+1) /2;
		for (j=0; j<=en; j++)
		{
			mpz_init (new_cache[idx+j]);
			new_ticky[idx+j] = 0;
		}
	}

	/* Now swap out the old and new */
	mpz_t* old_cache = c->cache;
	c->cache = new_cache;
	char* old_ticky = c->ticky;
	c->ticky = new_ticky;
	c->nmax = n;
	pthread_spin_unlock(&c->lock);

	if (old_cache) free(old_cache);
	if (old_ticky) free(old_ticky);
	return 0;
}

void i_one_d_cache_clear (i_cache *c)
{
	pthread_spin_lock(&c->lock);
	unsigned int i;
	for (i=0; i<c->nmax; i++)
	{
		c->ticky[i] = 0;
	}
	pthread_spin_unlock(&c->lock);
}

/* ======================================================================= */
/* Cache management */
/* pure cut-n-paste of he integer variant */

/** q_one_d_cache_check() -- check if mpq_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int q_one_d_cache_check (q_cache *c, unsigned int n)
{
	if ((n > c->nmax) || 0==c->nmax )
	{
		unsigned int newsize = 1.5*n+2;
		c->cache = (mpq_t *) realloc (c->cache, newsize * sizeof (mpq_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		unsigned int nstart = c->nmax+1;
		if (0 == c->nmax) nstart = 0;
		for (en = nstart; en < newsize; en++)
		{
			mpq_init (c->cache[en]);
			c->ticky[en] = 0;
		}
		c->nmax = newsize-1;
		return 0;
	}

	return (c->ticky[n]);
}

/* ======================================================================= */
/* Cache management */
/* Almost a cut-n-paste of above, but using fp instead */

/** fp_one_d_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int fp_one_d_cache_check (fp_cache *c, unsigned int n)
{
	pthread_spin_lock(&c->lock);
	if ((n <= c->nmax) && 0 != c->nmax)
	{
		int prec = c->precision[n];
		pthread_spin_unlock(&c->lock);
		return prec;
	}

	unsigned int newsize = 1.5*n+2;
	mpf_t* new_cache = (mpf_t *) malloc (newsize * sizeof(mpf_t));
	if (c->nmax) memcpy(new_cache, c->cache, (c->nmax+1) * sizeof(mpf_t));

	int* new_prec = (int *) malloc (newsize * sizeof(int));
	if (c->nmax) memcpy(new_prec, c->precision, (c->nmax+1) * sizeof(int));

	unsigned int en;
	unsigned int nstart = c->nmax+1;
	if (0 == c->nmax) nstart = 0;
	for (en = nstart; en < newsize; en++)
	{
		mpf_init (new_cache[en]);
		new_prec[en] = 0;
	}

	/* Now swap out the old and new */
	mpf_t* old_cache = c->cache;
	c->cache = new_cache;
	int* old_prec = c->precision;
	c->precision = new_prec;
	c->nmax = newsize-1;
	pthread_spin_unlock(&c->lock);

	if (old_cache) free(old_cache);
	if (old_prec) free(old_prec);
	return 0;
}

void fp_one_d_cache_clear (fp_cache *c)
{
	unsigned int i;
	pthread_spin_lock(&c->lock);
	for (i=0; i<c->nmax; i++)
	{
		c->precision[i] = 0;
	}
	pthread_spin_unlock(&c->lock);
}

/* ======================================================================= */
/** fp_triangle_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the
 *  value in the cache; else it returns zero.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int fp_triangle_cache_check (fp_cache *c, unsigned int n, unsigned int k)
{
	if ((n > c->nmax) || 0==c->nmax )
	{
		if (0 == n) n = 1;
		unsigned int newsize = (n+1)*(n+2)/2;
		c->cache = (mpf_t *) realloc (c->cache, newsize * sizeof (mpf_t));
		c->precision = (int *) realloc (c->precision, newsize * sizeof (int));

		unsigned int en;
      unsigned int nstart = c->nmax + 1;
      if (0 == c->nmax) nstart = 0;
		for (en = nstart; en <= n; en++)
		{
			unsigned int j;
			unsigned int idx = en * (en+1) /2 ;
			for (j=0; j<=en; j++)
			{
				mpf_init (c->cache[idx+j]);
				c->precision[idx+j]=0;
			}
		}
		c->nmax = n;
		return 0;
	}
	unsigned int idx = n * (n+1) /2 ;
	return c->precision[idx+k];
}

/* ======================================================================= */
/* Cache management */
/* A cut-n-paste of above, but using cpx instead */

/** cpx_one_d_cache_check() -- check if cpx_t value is in the cache
 *  If there is a cached value, this returns the precision of the
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int cpx_one_d_cache_check (cpx_cache *c, unsigned int n)
{
	pthread_spin_lock(&c->lock);
	if ((n <= c->nmax) && 0 != c->nmax)
	{
		int prec = c->precision[n];
		pthread_spin_unlock(&c->lock);
		return prec;
	}

	unsigned int newsize = 1.5*n+2;
	cpx_t* new_cache = (cpx_t *) malloc (newsize * sizeof(cpx_t));
	if (c->nmax) memcpy(new_cache, c->cache, (c->nmax+1) * sizeof(cpx_t));

	int* new_prec = (int *) malloc (newsize * sizeof(int));
	if (c->nmax) memcpy(new_prec, c->precision, (c->nmax+1) * sizeof(int));

	unsigned int en;
	unsigned int nstart = c->nmax+1;
	if (0 == c->nmax) nstart = 0;
	for (en = nstart; en < newsize; en++)
	{
		cpx_init (new_cache[en]);
		new_prec[en] = 0;
	}

	/* Now swap out the old and new */
	cpx_t* old_cache = c->cache;
	c->cache = new_cache;
	int* old_prec = c->precision;
	c->precision = new_prec;
	c->nmax = newsize-1;
	pthread_spin_unlock(&c->lock);

	if (old_cache) free(old_cache);
	if (old_prec) free(old_prec);
	return 0;
}

void cpx_one_d_cache_clear (cpx_cache *c)
{
	unsigned int i;
	pthread_spin_lock(&c->lock);
	for (i=0; i<c->nmax; i++)
	{
		c->precision[i] = 0;
	}
	pthread_spin_unlock(&c->lock);
}

/* =============================== END OF FILE =========================== */


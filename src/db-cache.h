/* 
 * db-cache.h
 *
 * File cache for pre-computed bignum values. Stores values,
 * with indicated precision, so we don't waste too much compute 
 * time.
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

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * fp_cache_put -- put mpf_t value into the database.
 * @prec: number of decimal places of precision to store.
 * @idx:  index under which to store the value
 */
void fp_cache_put (const char * dbname, const mpf_t val, int idx, int nprec);

/**
 * fp_cache_get -- get mpf_t from database
 * Returns 0 if no value in the database, or if the value in the
 * database has fewer than nprec digits. Thus, nprec is a minimum
 * requirement.
 */ 
int fp_cache_get (const char * dbname, mpf_t val, int idx, int nprec);

#ifdef  __cplusplus
};
#endif


/**
 * mp-zeroiso.c
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

#include <stdbool.h>
#include <stdlib.h>

#include "mp-complex.h"
#include "mp-zeroiso.h"

/* =============================================== */

// Bounding box
typedef struct box
{
	cpx_t boxll;
	cpx_t boxur;
	struct box* next;
} box_t;

static box_t* box_new(box_t* head, cpx_t boxll, cpx_t boxur)
{
	box_t* b = (box_t*) malloc(sizeof(box_t));
	cpx_init(b->boxll);
	cpx_init(b->boxur);
	cpx_set(b->boxll, boxll);
	cpx_set(b->boxur, boxur);
	b->next = head;
	head = b;
	return b;
}

static box_t* box_delete(box_t* head)
{
	cpx_clear(head->boxll);
	cpx_clear(head->boxur);
	box_t* next = head->next;
	free(head);
	return next;
}

// Split box into four
static box_t* box_split(box_t* head)
{
	cpx_t center;
	cpx_init(center);
	cpx_add(center, head->boxll, head->boxur);
	cpx_div_ui(center, center, 2);

	cpx_t cll;
	cpx_init(cll);
	cpx_t cur;
	cpx_init(cur);

	// The LL and UR quadrants.
	box_t* orig = head;
	cpx_set(cll, orig->boxll);  // ??? Compiler complains if I don't use indirection.
	head = box_new(head, cll, center);
	cpx_set(cur, orig->boxur);
	head = box_new(head, center, cur);

	// The UL quadrant
	mpf_set(cll[0].re, orig->boxll[0].re);
	mpf_set(cll[0].im, center[0].im);

	mpf_set(cur[0].re, center[0].re);
	mpf_set(cur[0].im, orig->boxur[0].im);

	head = box_new(head, cll, cur);

	// Clobber the original for the LR quadrant
	mpf_set(orig->boxll[0].re, center[0].re);
	// mpf_set(orig->boxll[0].im, unchanged.
	// mpf_set(orig->boxur[0].re, unchanged.
	mpf_set(orig->boxur[0].im, center[0].im);

	cpx_clear(center);
	cpx_clear(cll);
	cpx_clear(cur);
	return head;
}

static void box_midpoint(box_t* box, cpx_t center)
{
	cpx_add(center, box->boxll, box->boxur);
	cpx_div_ui(center, center, 2);
}

static void box_radius(box_t* box, mpf_t radius)
{
	mpf_t rim;
	mpf_init(rim);
	mpf_sub(radius, box->boxur[0].re, box->boxll[0].re);
	mpf_sub(rim, box->boxur[0].im, box->boxll[0].im);

	// The larger of the two.
	// mpf_cmp(a,b) is +1 if a>b
	if (0 < mpf_cmp(rim, radius))
		mpf_set(radius, rim);

	mpf_mul_ui(radius, radius, 3);
	mpf_div_ui(radius, radius, 4);
	mpf_clear(rim);
}

/* =============================================== */

// Test function. Implemented according to what the paper says.
// Offset should be zero or one.
static void test_fun(mpf_t bound,
              void (*poly)(cpx_t f, int deriv, cpx_t z, void* args),
              int degree, cpx_t center, mpf_t radius, int offset, void* args)
{
	cpx_t eval;
	cpx_init(eval);

	mpf_t term, rk, fact;
	mpf_init(term);
	mpf_init(rk);
	mpf_init(fact);

	mpf_set_ui(fact, 1);
	mpf_set_ui(rk, 1);

	mpf_set_ui(bound, 0);
	for (int k=1; k<=degree-offset; k++)
	{
		mpf_mul_ui(fact, fact, k);
		mpf_mul(rk, rk, radius);

		poly(eval, k+offset, center, args);
		cpx_abs(term, eval);
		mpf_mul(term, term, rk);
		mpf_div(term, term, fact);
		mpf_add(bound, bound, term);
	}

	// Divide by norm
	poly(eval, offset, center, args);
	cpx_abs(term, eval);
	mpf_div(bound, bound, term);

	cpx_clear(eval);
	mpf_clear(fact);
	mpf_clear(rk);
	mpf_clear(term);
}

/* =============================================== */

static mpf_t one;
static mpf_t hsqrt2;
static bool init = false;

// Test pedicate. Implemented according to what the paper says.
static bool test_predicate(
              void (*poly)(cpx_t f, int deriv, cpx_t z, void* args),
              int degree, cpx_t center, mpf_t radius, int offset, void* args)
{
	if (false == init)
	{
		mpf_init(one);
		mpf_set_ui(one, 1);
		mpf_init(hsqrt2);
		mpf_sqrt_ui(hsqrt2, 2);
		mpf_div_ui(hsqrt2, hsqrt2, 2);
		init = true;
	}
	mpf_t est;
	mpf_init(est);

	bool test = false;
	if (0 == offset)
	{
		test_fun(est, poly, degree, center, radius, 0, args);
		// mpf_cmp(a,b) is +1 if a>b
		test = (0 < mpf_cmp(one, est));
	}
	else
	{
		mpf_t rad;
		mpf_init(rad);
		mpf_mul_ui(rad, radius, 4*degree);

		test_fun(est, poly, degree, center, radius, 1, args);
		// mpf_cmp(a,b) is +1 if a>b
		test = (0 < mpf_cmp(hsqrt2, est));

		mpf_clear(rad);
	}

	mpf_clear(est);
	return test;
}

/* =============================================== */

// Return true if two dists overlap.
// Disk centers are ca, cb, disk radii are ra, rb
bool disk_intersect(cpx_t ca, mpf_t ra, cpx_t cb, mpf_t rb)
{
	cpx_t diff;
	cpx_init(diff);
	mpf_t dist;
	mpf_init(dist);

	cpx_sub(diff, ca, cb);
	cpx_abs(dist, diff);
	mpf_sub(dist, dist, ra);
	// mpf_cmp(a,b) is +1 if a>b
	bool rc = (0 < mpf_cmp(rb, dist));

	cpx_clear(diff);
	mpf_clear(dist);

	return rc;
}

/* =============================================== */

/**
 * Implements
 *    Michael Sagraloff, Chee K. Yap, "A Simple But Exact and Efficient
 *    Algorithm for Complex Root Isolation" (2011)
 *    https://cs.nyu.edu/exact/doc/complex.pdf
 */
int cpx_isolate_roots(
              void (*poly)(cpx_t f, int deriv, cpx_t z, void* args),
              int degree,
              cpx_t boxll, cpx_t boxur,
              cpx_t* centers, mpf_t* radii,
              void* args)
{
	cpx_t midpoint;
	cpx_init(midpoint);

	mpf_t radius;
	mpf_init(radius);

	int nfound = 0;

	// Avoid user error with swapped box coordinates
	if (0 < mpf_cmp(boxll[0].re, boxur[0].re))
	{
		mpf_set(radius, boxll[0].re);
		mpf_set(boxll[0].re, boxur[0].re);
		mpf_set(boxur[0].re, radius);
	}
	if (0 < mpf_cmp(boxll[0].im, boxur[0].im))
	{
		mpf_set(radius, boxll[0].im);
		mpf_set(boxll[0].im, boxur[0].im);
		mpf_set(boxur[0].im, radius);
	}

	box_t* head = box_new(NULL, boxll, boxur);
	while (NULL != head)
	{
		box_midpoint(head, midpoint);
		box_radius(head, radius);

		// First test; if it passes, the box is empty and discard it.
		bool result = test_predicate(poly, degree, midpoint, radius, 0, args);
		if (result)
		{
			head = box_delete(head);
			continue;
		}

		// Second test; if it fails, split box into four and try again.
		result = test_predicate(poly, degree, midpoint, radius, 1, args);
		if (false == result)
		{
			head = box_split(head);
			continue;
		}

		// Found a box with one root in it.
		mpf_mul_ui(radius, radius, 2*degree);

		bool replaced = false;
		for (int n=0; n<=nfound; n++)
		{
			if (disk_intersect(centers[n], radii[n], midpoint, radius))
			{
				// mpf_cmp(a,b) is +1 if a>b
				if (0 < mpf_cmp(radii[n], radius))
				{
					cpx_set(centers[n], midpoint);
					mpf_set(radii[n], radius);
				}
				replaced = true;
				break;
			}
		}
		if (false == replaced)
		{
			cpx_set(centers[nfound], midpoint);
			mpf_set(radii[nfound], radius);
			nfound++;
		}

		// We are done with this one.
		head = box_delete(head);
	}

	cpx_clear(midpoint);
	mpf_clear(radius);
	return nfound;
}

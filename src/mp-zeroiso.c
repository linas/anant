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

box_t* box_new(box_t* head, cpx_t boxll, cpx_t boxur)
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

box_t* box_delete(box_t* head)
{
	cpx_clear(head->boxll);
	cpx_clear(head->boxur);
	box_t* next = head->next;
	free(head);
	return next;
}

// Split box into four
box_t* box_split(box_t* head)
{
	cpx_t center;
	cpx_init(center);
	cpx_add(center, head->boxll, head->boxur);
	cpx_div_ui(center, center, 2);

	// The LL and UR quadrants.
	box_t* orig = head;
	head = box_new(head, orig->boxll, center);
	head = box_new(head, center, orig->boxur);

	// The UL quadrant
	cpx_t cll;
	cpx_init(cll);
	mpf_set(cll[0].re, orig->boxll[0].re);
	mpf_set(cll[0].im, center[0].im);

	cpx_t cur;
	cpx_init(cur);
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

void box_midpoint(box_t* box, cpx_t center)
{
	cpx_add(center, box->boxll, box->boxur);
	cpx_div_ui(center, center, 2);
}

void box_radius(box_t* box, mpf_t radius, mpf_t rim)
{
	mpf_sub(radius, box->boxll[0].re, box->boxur[0].re);
	mpf_sub(rim, box->boxll[0].im, box->boxur[0].im);

	// The larger of the two.
	if (0 < mpf_cmp(radius, rim))
		mpf_set(radius, rim);

	mpf_mul_ui(radius, radius, 3);
	mpf_div_ui(radius, radius, 4);
}

/**
 * Implements
 *    Michael Sagraloff, Chee K. Yap, "A Simple But Exact and Efficient
 *    Algorithm for Complex Root Isolation" (2011)
 *    https://cs.nyu.edu/exact/doc/complex.pdf
 */
int cpx_isolate_roots(
              void (*poly)(cpx_t f, int deriv, cpx_t z, void* args),
              int degree,
              cpx_t boxll, cpx_t boxrr,
              cpx_t* centers, mpf_t* radii,
              void* args)
{
	cpx_t midpoint;
	cpx_init(midpoint);

	mpf_t radius;
	mpf_init(radius);
	mpf_t tmp;
	mpf_init(tmp);

	int nfound = 0;

	box_t* head = box_new(NULL, boxll, boxrr);
	while (NULL != head)
	{
		box_midpoint(head, midpoint);
		box_radius(head, radius, tmp);
	}

	cpx_clear(midpoint);
	mpf_clear(radius);
	mpf_clear(tmp);
	return nfound;
}

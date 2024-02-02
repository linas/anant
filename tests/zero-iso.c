/*
 * zero-iso.c
 *
 * Debug unit test for zero isolation.
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

#include <stdio.h>
#include "mp-zeroiso.h"

/* ==================================================================== */

// Simple degree two polynomial, (x-1)(x-1/2)
void poly2(cpx_t f, int deriv, cpx_t z, void* args)
{
	cpx_t zn;
	cpx_init(zn);
	if (0 == deriv)
	{
		// + 1/2
		cpx_set_ui(f, 1, 0);
		cpx_div_ui(f, f, 2);

		// - 3x/2
		cpx_times_ui(zn, z, 3);
		cpx_div_ui(zn, zn, 2);
		cpx_sub(f, f, zn);

		// + x^2
		cpx_mul(zn, z, z);
		cpx_add(f, f, zn);
	}
	else if (1 == deriv)
	{
		// - 3/2
		cpx_set_ui(f, 3, 0);
		cpx_div_ui(f, f, 2);
		cpx_neg(f, f);

		// +2x
		cpx_times_ui(zn, z, 2);
		cpx_add(f, f, zn);
	}
	else if (2 == deriv)
	{
		cpx_set_ui(f, 2, 0);
	}
	else
	{
		fprintf(stderr, "Error: unexpected derivative %d\n", deriv);
		cpx_set_ui(f, 0, 0);
	}

	//printf("Poly call der= %d center= %f %f\n",
	//	deriv, cpx_get_re(z), cpx_get_im(z));
	cpx_clear(zn);
}

/* ==================================================================== */

// Degree n polynomial, (x-1)(x-1/2)(x-1/3)(x-1/4) etc.
// Computed recursively, so slow-ish for high orders.
void polyn(cpx_t f, int deriv, cpx_t z, void* args)
{
	int order = *((int*) args);

	cpx_t zn;
	cpx_init(zn);
	if (0 == deriv)
	{
		// Compute (z - 1/order)
		cpx_set_ui(zn, 1, 0);
		cpx_div_ui(zn, zn, order);
		cpx_sub(zn, z, zn);
		cpx_set_ui(f, 0, 1);

		// Recurse to get lower orders
		order --;
		if (0 < order)
			polyn(f, 0, z, &order);

		// Multiply (z - 1/order) * poly(lower orders)
		cpx_mul(f, f, zn);
	}
	else if (1 == order)
	{
		// Deriv of (z - 1/order) is just one, so we are done.
		cpx_set_ui(f, 0, 1);
	}
	else if (deriv <= order)
	{
		// Recurse to get lower orders
		// Compute (z - 1/order)
		cpx_set_ui(zn, 1, 0);
		cpx_div_ui(zn, zn, order);
		cpx_sub(zn, z, zn);
		cpx_set_ui(f, 0, 1);

		order --;
		polyn(f, deriv, z, &order);

		// Multiply (z - 1/order) * deriv(lower orders)
		cpx_mul(f, f, zn);

		// Deriv of (z - 1/order) is just one, so add lower
		polyn(zn, deriv-1, z, &order);
		cpx_add(f, f, zn);
	}
	else
	{
		fprintf(stderr, "Error:  derivative %d greater than order %d\n",
			deriv, order);
		cpx_set_ui(f, 0, 0);
	}

	cpx_clear(zn);
}

/* ==================================================================== */

int main (int argc, char * argv[])
{
	cpx_t boxll, boxur;
	cpx_init(boxll);
	cpx_init(boxur);

#define DEG 20
	cpx_t centers[DEG];
	mpf_t radii[DEG];
	for (int i=0; i<DEG; i++)
	{
		cpx_init(centers[i]);
		mpf_init(radii[i]);
	}

	cpx_set_ui(boxur, 2, 2);
	cpx_neg(boxll, boxur);

	// Simple polynomial with two roots.
	int nfound = cpx_isolate_roots(poly2, 2, boxll, boxur, centers, radii, NULL);

	printf("Found %d disks\n", nfound);
	for (int i=0; i<nfound; i++)
	{
		printf("Disk %d center= %f %f radius= %f\n", i,
			cpx_get_re(centers[i]),
			cpx_get_im(centers[i]),
			mpf_get_d(radii[i]));
	}

	// Degree n polynomial with roots colliding
	for (int degree=2; degree<15; degree++)
	{
		nfound = cpx_isolate_roots(polyn, degree, boxll, boxur,
			centers, radii, &degree);

		printf("Degree %d found %d disks\n", degree, nfound);
		for (int i=0; i<nfound; i++)
		{
			printf("Disk %d center= %f %f radius= %f\n", i,
				cpx_get_re(centers[i]),
				cpx_get_im(centers[i]),
				mpf_get_d(radii[i]));
		}
		printf("----\n");
	}

	return 0;
}

/* =============================== END OF FILE =========================== */

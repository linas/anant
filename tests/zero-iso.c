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
		fprintf(stderr, "Error: unexpected erivative %d\n", deriv);
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

	int nfound = cpx_isolate_roots(poly2, 2, boxll, boxur, centers, radii, NULL);

	printf("Found %d disks\n", nfound);
	for (int i=0; i<nfound; i++)
	{
		printf("Disk %d center= %f %f radius= %f\n", i,
			cpx_get_re(centers[i]),
			cpx_get_im(centers[i]),
			mpf_get_d(radii[i]));
	}

	return 0;
}

/* =============================== END OF FILE =========================== */

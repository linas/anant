/* 
 * db-merge.c
 *
 * Merge two file caches containing pre-computed bignum values. 
 * The value with the highest precision is kept.
 *
 * Linas Vepstas July 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#include "db-cache.h"

int
main (int argc, char * argv[])
{
	if (argc<4)
	{
		fprintf (stderr,"Usage: %s <db> <from> <to>\n", argv[0]);
		exit (1);
	}
	char * db = argv[1];
	int from = atoi (argv[2]);
	int maxidx = atoi (argv[3]);

	/* Compute number of binary bits this corresponds to. */
	int prec = 100;
	double v = ((double) prec) *log(10.0) / log(2.0);
	int bits = (int) (v + 100);
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);

	printf ("printout of %s up to n=%d\n", db, maxidx);

	mpf_t val;
	mpf_init (val);
	int n;
	int lastprec=-1;
	int lastn = 0;
	for (n=from; n<= maxidx; n++)
	{
		int prec = fp_cache_get (db, val, n, 10);
		if (prec != lastprec)
		{
			if (-1 != lastprec)
			{
				if (0 != lastn)
				{
					printf ("range %d : %d is to precision %d\n", lastn, n-1, lastprec);
				}
				else
				{
					printf ("range %d : %d is missing\n", lastn, n-1);
				}
			}
			lastprec = prec;
			lastn = n;
		}
	}
	if ((-1 != lastprec) && (0 != lastn))
	{
		printf ("range %d : %d is to precision %d\n", lastn, n-1, lastprec);
	}
	else
	{
		printf ("range %d : %d is missing\n", lastn, n-1);
	}
	return 0;
}

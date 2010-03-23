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
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "db-cache.h"
#include "mp-misc.h"

int
main (int argc, char * argv[])
{
	if (argc<5)
	{
		fprintf (stderr,"Usage: %s <out-db> <in-dba> <in-dbb> <from> <to>\n", argv[0]);
		exit (1);
	}
	char * dbout = argv[1];
	char * dbina = argv[2];
	char * dbinb = argv[3];
	int from = atoi (argv[4]);
	int too = atoi (argv[5]);

	/* Compute number of binary bits this corresponds to. */
	int prec = 10000;
	double v = ((double) prec) *log(10.0) / log(2.0);
	int bits = (int) (v + 100);
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);

	printf ("Merging %s and %s into %s from %d to %d\n", dbina, dbinb, dbout, from, too);

	mpf_t vala, valb;
	mpf_init (vala);
	mpf_init (valb);
	int n;
	for (n=from; n<= too; n++)
	{
		int preca = fp_cache_get (dbina, vala, n, 10);
		int precb = fp_cache_get (dbinb, valb, n, 10);

		if (0 >= preca && 0 >= precb) continue;
		
		if (precb < preca)
		{
			if (strcmp (dbina, dbout))
			{
				fp_cache_put (dbout, vala, n, preca);
				printf ("%d to %d from %s\t", n, preca, dbina);
				fp_prt ("", vala);
				printf ("\n");
			}
		}
		else
		{
			if (strcmp (dbinb, dbout))
			{
				fp_cache_put (dbout, valb, n, precb);
				printf ("%d to %d from %s\t", n, precb, dbinb);
				fp_prt ("", valb);
				printf ("\n");
			}
		}
	}
	return 0;
}

/*
 * mp-gkw.c
 *
 * FUNCTION:
 * Compute matrix elts of the GKW operator.
 *
 * HISTORY:
 * Linas Jan 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp-binomial.h"
#include "mp-misc.h"
#include "mp-zeta.h"


// Return the matrix element for the matrix element G_mp of the GKW
// operator, expanded at the x=1 location.
void
gkw(mpf_t acc, int m, int p, unsigned int prec)
{
	mpf_t one, term, fbin;
	mpz_t bin;
	int k;

	mpf_init (term);
	mpf_init (one);
	mpf_init (fbin);

	mpz_init (bin);

	// long double acc = 0.0L;
	mpf_set_ui (acc, 0);
	mpf_set_ui (one, 1);

	for (k=0; k<=p; k++)
	{
		// long double term = zetam1 (k+m+2);
		fp_zeta (term, k+m+2, prec);
		mpf_sub(term, term, one);

		// term *= binomial (m+k+1,m);
		i_binomial (bin, m+k+1, m);
		mpf_set_z (fbin, bin);
		mpf_mul (term, term, fbin);

		// term *= binomial (p,k);
		i_binomial (bin, p, k);
		mpf_set_z (fbin, bin);
		mpf_mul (term, term, fbin);

		if (k%2 == 0) mpf_add (acc, acc, term);
		else mpf_sub (acc, acc, term);
	}
}


// Return the continuous-valued version of the GKW operator.
// (the matrix elts occur at integer values)
// This implementation uses GMP multi-precision
void
gkw_smooth(mpf_t acc, double m, double p, unsigned int prec)
{
	mpf_t one, term, bin;
	int k;

	mpf_init (term);
	mpf_init (one);
	mpf_init (bin);

	cpx_t ess, zeta;
	cpx_init(ess);
	cpx_init(zeta);

	// long double acc = 0.0L;
	mpf_set_ui (acc, 0);
	mpf_set_ui (one, 1);

	int ip = (int) floor(p);
	for (k=0; k<= ip; k++)
	{
		double km2, km1;

printf ("duuude k=%d---------------------------------\n", k);
		// long double term = zetam1 (k+m+2);
		// fp_zeta (term, k+m+2, prec);
		km2 = k + m + 2.0L;
		cpx_set_d (ess, km2, 0.0);
		cpx_borwein_zeta(zeta, ess, prec);
		mpf_sub(term, zeta[0].re, one);
printf ("duuude k+m+2=%f\n", km2);
cpx_prt("cpx zeta = ", zeta);
printf("\n");
fp_prt("zeta -1 = ", term);
printf("\n");

		// term *= binomial (m+k+1,m);
		km1 = m+k+1.0L;
		fp_binomial_d (bin, km1, k+1);
		mpf_mul (term, term, bin);

		// term *= binomial (p,k);
		fp_binomial_d (bin, p, k);
		mpf_mul (term, term, bin);

		if (k%2 == 0) mpf_add (acc, acc, term);
		else mpf_sub (acc, acc, term);
	}
}

/* --------------------------- END OF LIFE ------------------------- */

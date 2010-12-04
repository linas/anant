
/** 
 * mp-polylog.c
 *
 * Implement Borwein-style polylogarithm.
 * Also implement the "periodic zeta" and 
 * the Hurwitz zeta function.
 * Also implements te Euler-Maclaurin expansion as well.
 *
 * As of 22 December 2006, seems to be fully functional
 * and correct, and passes tests. The range of convergence
 * is rather limited because of precision/rounding errors.
 * 
 * Copyright (C) 2006,2007 Linas Vepstas
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

#include <math.h>
#include <stdio.h>

#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-polylog.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* ============================================================= */
/**
 * polylog_borwein() -- Return polylog, Li_s(z) for estimator n.
 *
 * Uses the Borwein algorithm to estimate the polylog value.
 * A polynomial of order 2n is used to perform the estimation.
 * The value of n must be suitably large.
 * Muse have |z^2/(z-1)| < 3.7 in order for the thing to converge
 * to an accurate answer. 
 *
 * Appears to work well. Suggest n=31 for most cases,
 * should return answers accurate to 1e-16
 */
static void polylog_borwein (cpx_t plog, const cpx_t ess, const cpx_t zee, int norder, int prec)
{
	DECLARE_CPX_CACHE (bin_sum);
	mpz_t ibin;
	cpx_t s, z, ska, pz, acc, term, ck, bins;
	int k;

	mpz_init (ibin);
	cpx_init (s);
	cpx_init (z);
	cpx_init (ska);
	cpx_init (pz);
	cpx_init (acc);
	cpx_init (term);
	cpx_init (ck);
	cpx_init (bins);

	/* s = -ess */
	cpx_neg (s, ess);
	cpx_set (z, zee);

	/* first binomial summation term is 1 */
	cpx_set_ui (bins, 1, 0);
	cpx_one_d_cache_check (&bin_sum, 0);
	cpx_one_d_cache_store (&bin_sum, bins, 0, prec);

	/* ska = [1/(z-1)]^n */
	cpx_set (ska, z);
	cpx_sub_ui (ska, ska, 1, 0);
	cpx_recip (ska, ska);
	cpx_pow_ui (ska, ska, norder);
	
	cpx_set_ui (pz, 1, 0);
	cpx_set_ui (acc, 0, 0);
	cpx_set_ui (plog, 0, 0);

	for (k=1; k<=norder; k++)
	{
		cpx_mul(pz, pz, z);

		/* The inverse integer power */
		cpx_ui_pow_cache (term, k, s, prec);

		/* Put it together */
		cpx_mul (term, term, pz);
		cpx_add (acc, acc, term);

		/* Compute the binomial sum */
		i_binomial (ibin, norder, k);
		mpf_set_z (term[0].re, ibin);
		mpf_set_ui (term[0].im, 0);
		cpx_mul (term, term, pz);
		
		if (k%2)
		{
			cpx_sub (bins, bins, term);
		}
		else
		{
			cpx_add (bins, bins, term);
		}

		/* Stow the binomial sum away in an array;
		 * we'll need to reference this in reverse order later.
		 */
		cpx_one_d_cache_check (&bin_sum, k);
		cpx_one_d_cache_store (&bin_sum, bins, k, prec);
	}

	for (k=norder+1; k<=2*norder; k++)
	{
		cpx_mul(pz, pz, z);

		/* The inverse integer power */
		cpx_ui_pow_cache (term, k, s, prec);
		cpx_mul (term, term, pz);

		/* Fetch binomial sum from the array */
		cpx_one_d_cache_fetch (&bin_sum, bins, 2*norder-k);
		cpx_mul (term, term, bins);

		/* Put it together */
		cpx_add (plog, plog, term);
	}

	cpx_mul (plog, plog, ska);
	if (norder%2)
	{
		cpx_sub (plog, acc, plog);
	}
	else
	{
		cpx_add (plog, acc, plog);
	}
	
	cpx_clear (s);
	cpx_clear (z);
	cpx_clear (ska);
	cpx_clear (pz);
	cpx_clear (acc);
	cpx_clear (term);
	cpx_clear (ck);
	cpx_clear (bins);
	mpz_clear (ibin);

	cpx_one_d_cache_clear(&bin_sum);
}

/* ============================================================= */

/* polylog_get_zone -- return | z^2 / (z-1) |^2
 *
 * The value of | z^2 / (z-1) | is used to determine
 * the convergence zone.
 */
inline static double polylog_get_zone (double zre, double zim)
{
	double den = 1.0 / ((zre-1.0)*(zre-1.0) + zim*zim);
	double sre = zre*zre - zim*zim;
	double sim = 2.0*zre*zim;
	double fre = sre * (zre-1.0) + zim*sim;
	double fim = sim * (zre-1.0) - zim*sre;
	den = (fre*fre + fim*fim)*den*den;

	return den;
}

inline static double polylog_modsq (const cpx_t zee)
{
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	

	double den = zre*zre + zim*zim;

	return den;
}

/*
 * polylog_terms_est() -- estimate number of terms needed 
 * in the polylog summation in order to keep the error
 * to be less than 10^-prec.
 */
static int polylog_terms_est (const cpx_t ess, const cpx_t zee, int prec)
{
	double fterms = 2.302585 * prec;  /* log(10) */

	/* Estimate for the gamma. A slightly better estimate
	 * can be obtains for sre negative but still small. 
	 */
	double gamterms;
	double sre = mpf_get_d (ess[0].re);	
	double sim = mpf_get_d (ess[0].im);	
	if (0.0 > sim) sim = -sim;
	if (0.0 < sre) {
		gamterms = 0.5*M_PI*sim;
	} else {
		gamterms = M_PI*sim;
	}
	/* XXX TODO replace lgamma with stirling approx slog(s)-s for better
	 * performance and speed */
	gamterms -= lgamma(sre); 

	/*
	 * If lngamma is divergent, then sre is a negative integer,
	 * and the approximation is in fact an exact result. Thus,
	 * we need only set nterms to be minus the real part of s.
	 *
	 * XXXX this is wrong, this would hold only for
	 * sim = 0. But if someone passes in sre=-n and sim != 0,
	 * then this does the wrong thing :-( Boooo...
	 */
	if ((-10123>gamterms) || (10123 < gamterms))
	{
		return (int) (-sre+3.0);
	}
	
	fterms += gamterms;

	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	
	double cterms = 0.0;
	if (0.0 < zre)
	{
		double mod = zre*zre + zim*zim;
		if (1.0 >mod)
		{
			double mod = (zre-1.0)*(zre-1.0) + zim*zim;
			cterms = -0.5 * log (mod);
		}
		else
		{
			cterms = 0.5 * log (mod);
			if (0.0 > zim) zim = -zim;
			cterms -= log (zim);
		}
		fterms += cterms;
	}

	/* den = | z^2/(z-1)|^2 */
	double den = polylog_get_zone (zre, zim);

	/* fterms may become negative -- a negative value means 
	 * it will never converge. Caller must test for negative value.
	 */
	fterms /= -0.5*log(den) + 1.386294361;  /* log 4 */
	int nterms = (int) (fterms+1.0);

#if 0
gamterms /=  -0.5*log(den) + 1.386294361;
cterms /=  -0.5*log(den) + 1.386294361;
printf ("# duude z= %g +i %g den=%g  prec=%d deno=%g nterms = %d gam=%g ct=%g\n", 
zre, zim, sqrt(den), prec,  -0.5*log(den) + 1.386294361, nterms, gamterms, cterms);
#endif

	return nterms;
}

/* ============================================================= */

static int recurse_away_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth);
  
static inline int polylog_recurse_duple (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int rc;
	cpx_t zsq, s, pp, pn;
	cpx_init (zsq);
	cpx_init (s);
	cpx_init (pp);
	cpx_init (pn);

	cpx_mul (zsq, zee, zee);
	cpx_set (s, ess);

#if 0
cpx_prt ("dupl zee= ", zee);
printf ("\n");
cpx_prt ("zsq= ", zsq);
printf ("\n");
#endif
	rc = recurse_away_polylog (pp, s, zsq, prec, depth);
	if (rc) goto bailout;

	cpx_neg (zsq, zee);
	rc = recurse_away_polylog (pn, s, zsq, prec, depth);
	if (rc) goto bailout;

	/* now, compute 2^{1-s} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_ui_pow (s, 2, s, prec);
	cpx_mul (plog, pp, s);

	cpx_sub (plog, plog, pn);
	
bailout:
	cpx_clear (s);
	cpx_clear (pp);
	cpx_clear (pn);
	cpx_clear (zsq);
	return rc;
}

static inline int polylog_recurse_triple (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int rc;
	cpx_t zcu, s, tr, pp, pu, pd;
	cpx_init (zcu);
	cpx_init (s);
	cpx_init (tr);
	cpx_init (pp);
	cpx_init (pu);
	cpx_init (pd);

	cpx_set (s, ess);

	cpx_mul (zcu, zee, zee);
	cpx_mul (zcu, zcu, zee);
	rc = recurse_away_polylog (pp, s, zcu, prec, depth);
	if (rc) goto bailout;

	/* tr = exp (i 2pi/3) = -1/2  + i sqrt(3)/2 */
	mpf_set_ui (tr[0].re, 1);
	mpf_div_ui (tr[0].re, tr[0].re, 2);
	mpf_neg (tr[0].re, tr[0].re);
	fp_half_sqrt_three (tr[0].im);
	
	cpx_mul (zcu, tr, zee);
	rc = recurse_away_polylog (pu, s, zcu, prec, depth);
	if (rc) goto bailout;

	cpx_mul (zcu, tr, zcu);
	rc = recurse_away_polylog (pd, s, zcu, prec, depth);
	if (rc) goto bailout;

	/* now, compute 3^{1-s} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_ui_pow (s, 3, s, prec);
	cpx_mul (plog, pp, s);

	cpx_sub (plog, plog, pu);
	cpx_sub (plog, plog, pd);
	
bailout:
	cpx_clear (s);
	cpx_clear (tr);
	cpx_clear (pp);
	cpx_clear (pu);
	cpx_clear (pd);
	cpx_clear (zcu);
	return rc;
}

/**
 * recurse_away_polylog() -- use duplication formula to extend domain
 *
 * Evaluate the polylog directly, if possible; else use the 
 * duplication formula to get into a region where its directly 
 * evaluable. The duplication formula is used to move away from
 * z=1, and the Hurwitz series at z=1 is *not* used.
 * 
 * Return a non-zero value if no value was computed.
 */
static int recurse_away_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int rc;
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	
	double mod = zre*zre + zim*zim;

	/* The algo will never converge when modulus >= 5 or so */
	if (25 < mod) return 1;

	/*
	 * Limit the depth of recursion to avoid run-away. Now
	 * that the algo is working well, this seems to almost
	 * never be needed (!?).
	 */
	if (9 < depth)
	{
		fprintf (stderr, "excessive recursion (away) at z=%g+ i%g\n", zre, zim);
		return 1;
	}
	depth ++;

	/*
	 * The Borwein algo seems to always be faster than direct 
	 * summation, even when the direct-sum region is made quite 
	 * small, e.g. even when it is of radius less than 1/4.
	 * Never use direct summation.
	 */
#if 0
	if (0.0625>mod)
	{
		cpx_polylog_sum (plog, ess, zee, prec);
		return 0;
	}
#endif

	/* The zone of convergence for the Borwein algorithm is 
	 * |z^2/(z-1)| < 3.  If z is within this zone, then all is 
	 * well. If not, use the duplication formula to make 
	 * recursive calls, until the leaves of the recursion 
	 * are in this zone. 
	 *
	 * The algo seems to be more precise (!??) and have less 
	 * trouble when an even smaller bound is used, e.g.
	 * when |z^2/(z-1)| < 1.25. The non-recursive algo seems
	 * to choke up when it gets too close to z=1.
	 */

	/* den = | z^2/(z-1)|^2 */
	double den = polylog_get_zone (zre, zim);

	/* nterms = number of terms needed for computation */
	int nterms = polylog_terms_est (ess, zee, prec);

	/* To carry out the computation, an internal precision is needed 
	 * that is a bit higher than what the user asked for. This internal
	 * precision depends on the degree of the approximating polynomial.
	 * Thus, look at the available bits of precision, and decide if
	 * the calculation can be performed in that way. 
	 *
	 * The degree of internal precision available limits the largest
	 * effective order of the apprximating polynomial that can be used.
	 * Its pointless/erroneous to try to use a polynomial of degree
	 * more than "maxterms".
	 */
	int nbits = mpf_get_default_prec();
	int maxterms = nbits - (int) (3.321928095 *prec); /* log(10) / log(2) */

	// printf ("invoke-away, z=%g +i %g  den=%g nterms=%d, maxterms=%d\n", zre, zim, den, nterms, maxterms);
	/* if (4> nterms) (i.e. nterms is negative), then the thing will
	 * never converge, and so subdivision is the only option.
	 * Uhh, this should be equivalent to den>15, and we already subdivide
	 * for large den.
	 * 
	 * The algo seems to have some trouble near z=1 when 
	 * if (den>4) is used to decide subdivision. 
	 */
	if ((den > 1.5) || (maxterms < nterms))
	{
		// printf ("splitsville-away, z=%g +i %g  den=%g nterms=%d\n", zre, zim, den, nterms);
		rc = polylog_recurse_duple (plog, ess, zee, prec, depth);
		/* 
		 * The angle-tripling recursion equation is not as effective 
		 * as the angle-doubling equation in pulling points into the
		 * zone of convergence. So we don't use it.
		 *
		 * rc = polylog_recurse_triple (plog, ess, zee, prec, depth);
		 */
		return rc;

		/* 
		 * Under no circumstances does it ever seem to work to bring
		 * distant points closer in by using the sqrt relation to 
		 * pull them in. The problem seems to be that distant points
		 * get pulled in close to z>=1, where they can't be evaluated 
		 * anyway, so this is no particular help.
		 * 
		 * rc = polylog_recurse_sqrt (plog, ess, zee, prec, depth);
		 */
	}
	
	/* Use the larger, adjusted internal precision discussed above
	 * in the final calculation.
	 */
	prec += (int) (0.301029996 * nterms) +1;
	polylog_borwein (plog, ess, zee, nterms, prec);
	return 0;
}

int cpx_polylog_away (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	int rc = recurse_away_polylog (plog, ess, zee, prec, 0);
	if (rc)
	{
		cpx_set_ui (plog, 0,0);
	}
	return rc;
}

/* ============================================================= */

static int recurse_towards_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth);

static inline int 
polylog_recurse_sqrt (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int rc;
	cpx_t zroot, s, pp, pn;
	cpx_init (zroot);
	cpx_init (s);
	cpx_init (pp);
	cpx_init (pn);

	cpx_sqrt (zroot, zee, prec);
	cpx_set (s, ess);

#if 0
cpx_prt ("zee= ", zee);
printf ("\n");
cpx_prt ("zroot= ", zroot);
printf ("\n");
#endif
	rc = recurse_towards_polylog (pp, s, zroot, prec, depth);
	if (rc) goto bailout;

	cpx_neg (zroot, zroot);
	rc = recurse_towards_polylog (pn, s, zroot, prec, depth);
	if (rc) goto bailout;

	cpx_add (plog, pp, pn);

	/* now, compute 2^{s-1} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_ui_pow (s, 2, s, prec);
	cpx_mul (plog, plog, s);

bailout:
	cpx_clear (s);
	cpx_clear (pp);
	cpx_clear (pn);
	cpx_clear (zroot);
	return rc;
}

/*
 * polylog_invert -- implement the polylog inversion formula
 * Implement the following inversion formula for polylog:
 * (1-e^{2pi is}) Li_s(z) =  e^{i\pi s} (2pi i)^s / Gamma(s) 
 *           (zeta(1-s, ln z/(2pi i) -e^{ipi s} zeta(1-s, 1- ln z/(2pi i)) 
 *
 * This formula appears to work well for both positive and negative 
 * half s-plane.
 */
static int 
polylog_invert_works(cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_t s, oz, tmp, ph, term, logz;
	cpx_init (s);
	cpx_init (oz);
	cpx_init (tmp);
	cpx_init (term);
	cpx_init (ph);
	cpx_init (logz);
	cpx_set (s, ess);
	cpx_recip (oz, zee);

	/* compute ph = e^{i pi s / 2} = i^s */
	cpx_times_mpf (tmp, s, twopi);
	cpx_div_ui (tmp,tmp, 4);
	cpx_times_i (tmp, tmp);
	cpx_exp (ph, tmp, prec);

	/* compute ln z/(2pi i) */
	cpx_log (logz, oz, prec);
	cpx_div_mpf (logz, logz, twopi);
	cpx_times_i (logz, logz);

	/* Place branch cut so that it extends to the right from z=1 */
	if (mpf_sgn(logz[0].re) < 0)
	{
		mpf_add_ui (logz[0].re, logz[0].re, 1);
	}

	/* zeta (1-s, ln z/(2pi i)) */
	cpx_ui_sub (tmp, 1, 0, s);
	// cpx_hurwitz_taylor (term, tmp, logz, prec);
	cpx_hurwitz_euler (term, tmp, logz, prec);
	
	/* plus e^{ipi s} zeta (1-s, 1-ln z/(2pi i)) */
	cpx_neg (logz, logz);
	cpx_add_ui (logz, logz, 1, 0);
	// cpx_hurwitz_taylor (tmp, tmp, logz, prec);
	cpx_hurwitz_euler (tmp, tmp, logz, prec);
	cpx_mul (tmp, tmp, ph);
	cpx_mul (tmp, tmp, ph);
	cpx_sub (term, term, tmp);

	/* (2pi)^s i^s zeta /gamma (s) */
	cpx_mul (term, term, ph);
	cpx_mpf_pow (tmp, twopi, s, prec);
	cpx_mul (term, term, tmp);
	cpx_gamma_cache (tmp, s, prec);
	cpx_div (term, term, tmp);

	cpx_set (plog, term);

	/* divide by (1-e^{2pi s}) */
	cpx_mul(ph, ph, ph);
	cpx_mul(ph, ph, ph);

	cpx_neg(ph,ph);
	cpx_add_ui (ph, ph, 1, 0);
	cpx_div (plog, plog, ph);

	cpx_clear (s);
	cpx_clear (oz);
	cpx_clear (logz);
	cpx_clear (tmp);
	cpx_clear (term);
	cpx_clear (ph);
	mpf_clear (twopi);
	return 0;
}

/* The following is an alternate version of the invert routine,
 * it works, too, for the upper half or the lower half plane.
 * The only difference between this and the above is that the 
 * gamma function reflection formula was used to put the gamma
 * on the top and not the bottom.
 */
static int 
polylog_invert(cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int redo = 0;

	/* create a cache of commonly re-used values */
	static int cache_prec = -1;
	static mpf_t twopi, otp, log_twopi;
	static cpx_t phase, scale, cache_ess, s;
	if (-1 == cache_prec)
	{
		mpf_init (twopi);
		mpf_init (otp);
		mpf_init (log_twopi);

		cpx_init (phase);
		cpx_init (scale);
		cpx_init (s);
		cpx_init (cache_ess);
		cpx_set_ui(cache_ess, 123123123, 321321321);
	}

	if (cache_prec != prec)
	{
		redo = 1;
		cache_prec = prec;
		mpf_set_prec (twopi, 3.322*prec +50);
		mpf_set_prec (otp, 3.322*prec +50);
		mpf_set_prec (log_twopi, 3.322*prec +50);

		cpx_set_prec (phase, 3.322*prec +50);
		cpx_set_prec (scale, 3.322*prec +50);
		cpx_set_prec (s, 3.322*prec +50);
		cpx_set_prec (cache_ess, 3.322*prec +50);

		fp_two_pi (twopi, prec);

		/* otp = -1/2pi */
		mpf_set_ui (otp, 1);
		mpf_neg (otp, otp);
		mpf_div (otp, otp, twopi);

		fp_log (log_twopi, twopi, prec);
	}

	cpx_t oz, tmp, logz;
	cpx_init (oz);
	cpx_init (tmp);
	cpx_init (logz);

	/* Recompute these values only if s differs from last time. */
	if(redo || !cpx_eq (ess, cache_ess, prec*3.322))
	{
		cpx_set (cache_ess, ess);
		cpx_ui_sub (s, 1, 0, ess);

		/* compute ph = e^{-i pi s / 2} = (-i)^s */
		cpx_times_mpf (tmp, s, twopi);
		cpx_div_ui (tmp, tmp, 4);
		cpx_times_i (tmp, tmp);
		cpx_neg (tmp, tmp);
		cpx_exp (oz, tmp, prec);
		cpx_mul (phase, oz, oz);

		/* gamma(s) (i)^s / (2pi)^s */
		cpx_gamma_cache (scale, s, prec);
		cpx_times_mpf (tmp, s, log_twopi);
		cpx_neg (tmp, tmp);
		cpx_exp (tmp, tmp, prec);
		cpx_mul (scale, scale, tmp);
		cpx_div (scale, scale, oz);
	}

	/* compute ln z/(2pi i) */
	cpx_set (oz, zee);
	cpx_log (logz, oz, prec);
	cpx_times_mpf (logz, logz, otp);
	cpx_times_i (logz, logz);

	/* Place branch cut so that it extends to the right from z=1 */
	if (mpf_sgn(logz[0].re) < 0)
	{
		mpf_add_ui (logz[0].re, logz[0].re, 1);
	}

	/* zeta (s, ln z/(2pi i)) */
	// cpx_hurwitz_taylor (plog, s, logz, prec);
	cpx_hurwitz_euler (plog, s, logz, prec);
	
	/* plus e^{-ipi s} zeta (s, 1-ln z/(2pi i)) */
	cpx_ui_sub (logz, 1, 0, logz);
	// cpx_hurwitz_taylor (tmp, s, logz, prec);
	cpx_hurwitz_euler (tmp, s, logz, prec);
	cpx_mul (tmp, tmp, phase);
	cpx_add (plog, plog, tmp);

	cpx_mul (plog, plog, scale);

	cpx_clear (oz);
	cpx_clear (logz);
	cpx_clear (tmp);
	return 0;
}

// #define NON_WORKING_INVERSION_ROUTINES
#ifdef NON_WORKING_INVERSION_ROUTINES
/* Implement the following inversion formula for polylog:
 * Li_s(z) = - e^{i\pi s} Li_s(1/z) 
 *           + (2pi i)^s zeta(1-s, ln z/(2pi i)) / Gamma (s)
 *
 * This polylog inversion formula "should" work well.
 * and it does, for the upper half s-plane. However, for the
 * lower-half s-plane, its broken, and the reason for this
 * brokenness is confusing, since its theoretically the same
 * formula as the one that works. Not clear what the gig is.
 */
static int 
polylog_invert_broken_for_lower_half_plane(cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_t s, oz, tmp, ph, term, logz;
	cpx_init (s);
	cpx_init (oz);
	cpx_init (tmp);
	cpx_init (term);
	cpx_init (ph);
	cpx_init (logz);
	cpx_set (s, ess);
	cpx_recip (oz, zee);

	/* compute ph = e^{i pi s / 2} = i^s */
	cpx_times_mpf (tmp, s, twopi);
	cpx_div_ui (tmp,tmp, 4);
	cpx_times_i (tmp, tmp);
	cpx_exp (ph, tmp, prec);

	/* - e^i\pi s Li_s(1/z) */
	int rc = recurse_towards_polylog (plog, s, oz, prec, depth);
	if (rc) goto bail;
	
	cpx_mul (plog, plog, ph);
	cpx_mul (plog, plog, ph);
	cpx_neg (plog, plog);

	/* compute ln z/(2pi i) */
	cpx_log (logz, oz, prec);
	cpx_div_mpf (logz, logz, twopi);
	cpx_times_i (logz, logz);

	/* Place branch cut so that it extends to the right from z=1 */
	if (mpf_sgn(logz[0].re) < 0)
	{
		mpf_add_ui (logz[0].re, logz[0].re, 1);
	}

	/* zeta (1-s, ln z/(2pi i)) */
	cpx_ui_sub (tmp, 1, 0, s);
	cpx_hurwitz_taylor (term, tmp, logz, prec);
	
	/* (2pi)^s i^s zeta /gamma (s) */
	cpx_mul (term, term, ph);

	cpx_mpf_pow (tmp, twopi, s, prec);
	cpx_mul (term, term, tmp);
	cpx_gamma_cache (tmp, s, prec);
	cpx_div (term, term, tmp);

	cpx_add (plog, plog, term);

bail:
	cpx_clear (s);
	cpx_clear (oz);
	cpx_clear (logz);
	cpx_clear (tmp);
	cpx_clear (term);
	cpx_clear (ph);
	mpf_clear (twopi);
	return rc;
}
#endif /* NON_WORKING_INVERSION_ROUTINES */

/*
 * cpx_polylog_sheet -- move to sheet N of polylog
 * The Nth sheet of Li_s(z) is given by
 *      (2pi i)^s zeta(1-s, ln z/(2pi i)) / Gamma (s)
 */
void 
cpx_polylog_sheet(cpx_t delta, const cpx_t ess, const cpx_t zee, int z0_dromy, int z1_dromy, int prec)
{
	if (0 == z1_dromy)
	{
		cpx_set_ui (delta, 0,0);
		return;
	}

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_t s, tmp, ph, q, norm;
	cpx_init (s);
	cpx_init (tmp);
	cpx_init (ph);
	cpx_init (q);
	cpx_init (norm);
	cpx_set (s, ess);
	cpx_set_ui (norm, 1, 0);

	/* Do the z0_dromy */
	if (z0_dromy) 
	{
		cpx_times_mpf (tmp, s, twopi);
		cpx_times_i (tmp, tmp);
		cpx_neg (tmp, tmp);
		if (z0_dromy > 0)
		{
			cpx_times_ui (tmp, tmp, z0_dromy);
		}
		else
		{
			cpx_times_ui (tmp, tmp, -z0_dromy);
			cpx_neg (tmp,tmp);
		}
		cpx_exp (norm, tmp, prec);
		if (z0_dromy%2) 
		{
			cpx_neg (norm, norm);
		}
	}

	/* Compute q = ln z/(2pi i) */
	cpx_log (q, zee, prec);
	cpx_div_mpf (q, q, twopi);
	cpx_times_i (q, q);
	cpx_neg (q,q);
	
	/* Place branch cut of the polylog so that it extends to the 
	 * right from z=1. This is the same as adding 2pi i to the value 
	 * of the log, if the value is in the lower half plane, so that
	 * we move the cut of the log to lie along the positive, not 
	 * negative axis. */
	if (mpf_sgn(q[0].re) < 0)
	{
		mpf_add_ui (q[0].re, q[0].re, 1);
	}

	/* Move to the n'th sheet; sheets of the log and the polylog
	 * are now one and the same thing. */
	if (0 < z1_dromy)
	{
		mpf_add_ui (q[0].re, q[0].re, z1_dromy);
	}
	else
	{
		cpx_neg (q, q);
		mpf_add_ui (q[0].re, q[0].re, -z1_dromy);

		/* and one more, for the loop */
		mpf_add_ui (q[0].re, q[0].re, 1);
	}

	/* Compute sum over 1/q^s = (ln z/(2pi i))^{s-1} */
	cpx_sub_ui (s, s, 1, 0);
	cpx_set_ui (delta, 0, 0);

	while (mpf_cmp_ui (q[0].re, 1) > 0)
	{
		mpf_sub_ui (q[0].re, q[0].re, 1);
		cpx_pow (tmp, q, s, prec);
		cpx_add (delta, delta, tmp);
	}

	cpx_add_ui (s, s, 1, 0);

	/* compute normalization */
	// XXX this code can be more optimized than this
	// current form is easier to validate vs. theory
	/* Compute ph = e^{i pi s / 2} = i^s */
	cpx_times_mpf (tmp, s, twopi);
	cpx_div_ui (tmp, tmp, 4);
	cpx_times_i (tmp, tmp);

	if (0> z1_dromy)
	{
		cpx_neg (tmp, tmp);
	}
	cpx_exp (ph, tmp, prec);

	/* (2pi)^s i^s (sum) /gamma (s) */
	cpx_mul (delta, delta, ph);
	cpx_mpf_pow (tmp, twopi, s, prec);
	cpx_mul (delta, delta, tmp);
	cpx_gamma_cache (tmp, s, prec);
	cpx_div (delta, delta, tmp);

	cpx_mul (delta, delta, norm);

	cpx_clear (s);
	cpx_clear (q);
	cpx_clear (tmp);
	cpx_clear (ph);
	cpx_clear (norm);
	mpf_clear (twopi);
}

void 
cpx_polylog_sheet_g0_action(cpx_t ph, const cpx_t ess, int direction, int prec)
{
	if (0 == direction)
	{
		cpx_set_ui (ph, 0,0);
		return;
	}

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_times_mpf (ph, ess, twopi);
	cpx_times_i (ph, ph);
	cpx_neg (ph, ph);
	if (direction > 0)
	{
		cpx_times_ui (ph, ph, direction);
	}
	else
	{
		cpx_times_ui (ph, ph, -direction);
		cpx_neg (ph, ph);
	}
	cpx_exp (ph, ph, prec);
	if (direction%2) 
	{
		cpx_neg (ph, ph);
	}
	mpf_clear (twopi);
}

void 
cpx_polylog_sheet_g1_action(cpx_t delta, const cpx_t ess, const cpx_t zee, int sheet, int direction, int prec)
{
	if (0 == direction)
	{
		cpx_set_ui (delta, 0,0);
		return;
	}

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_t s, tmp, ph, q;
	cpx_init (s);
	cpx_init (tmp);
	cpx_init (ph);
	cpx_init (q);
	cpx_set (s, ess);

	/* Compute q = ln z/(2pi i) */
	cpx_log (q, zee, prec);
	cpx_div_mpf (q, q, twopi);
	cpx_times_i (q, q);
	cpx_neg (q,q);
	
	/* Place branch cut of the polylog so that it extends to the 
	 * right from z=1. This is the same as adding 2pi i to the value 
	 * of the log, if the value is in the lower half plane, so that
	 * we move the cut of the log to lie along the positive, not 
	 * negative axis. */
	if (mpf_sgn(q[0].re) < 0)
	{
		mpf_add_ui (q[0].re, q[0].re, 1);
	}

	/* Move to the n'th sheet; sheets of the log and the polylog
	 * are now one and the same thing. */
	// XXX this is wrong, if zero is crossed. 
	// That is, this works correctly only if direction is +1 or -1 
	// or if z1_dromy is same sign as sheet.  But for now, this
	// restriction is enough. 
	int z1_dromy = sheet + direction;
	if (0 < z1_dromy)
	{
		mpf_add_ui (q[0].re, q[0].re, z1_dromy);
	}
	else
	{
		cpx_neg (q, q);
		mpf_add_ui (q[0].re, q[0].re, -z1_dromy);

		/* and one more, for the loop */
		mpf_add_ui (q[0].re, q[0].re, 1);
	}

	/* Compute sum over 1/q^s = (ln z/(2pi i))^{s-1} */
	cpx_sub_ui (s, s, 1, 0);
	cpx_set_ui (delta, 0, 0);

	if (0 > direction) direction = -direction;
	while (direction > 0)
	{
		mpf_sub_ui (q[0].re, q[0].re, 1);
		cpx_pow (tmp, q, s, prec);
		cpx_add (delta, delta, tmp);
		direction --;
	}

	cpx_add_ui (s, s, 1, 0);

	/* compute normalization */
	// XXX this code can be more optimized than this
	// current form is easier to validate vs. theory
	/* Compute ph = e^{i pi s / 2} = i^s */
	cpx_times_mpf (tmp, s, twopi);
	cpx_div_ui (tmp, tmp, 4);
	cpx_times_i (tmp, tmp);

	if (0> z1_dromy)
	{
		cpx_neg (tmp, tmp);
	}
	cpx_exp (ph, tmp, prec);

	/* (2pi)^s i^s (sum) /gamma (s) */
	cpx_mul (delta, delta, ph);
	cpx_mpf_pow (tmp, twopi, s, prec);
	cpx_mul (delta, delta, tmp);
	cpx_gamma_cache (tmp, s, prec);
	cpx_div (delta, delta, tmp);

	cpx_clear (s);
	cpx_clear (q);
	cpx_clear (tmp);
	cpx_clear (ph);
	mpf_clear (twopi);
}

/**
 * recurse_towards_polylog() -- use duplication formula to extend domain
 *
 * Evaluate the polylog directly, if possible; else use the 
 * duplication formula to get into a region where its directly 
 * evaluable. The duplication formula is used to move towards z=1
 * which is where the Hurwitz series at z=1 can be employed.
 * 
 * Return a non-zero value if no value was computed.
 */
static int recurse_towards_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	int rc;
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	
	double mod = zre*zre + zim*zim;

	/*
	 * Limit the depth of recursion to avoid run-away. Now
	 * that the algo is working well, this seems to almost
	 * never be needed (!?).
	 */
	if (5 < depth)
	{
		fprintf (stderr, "excessive recursion (to) at z=%g+ i%g\n", zre, zim);
		return 1;
	}
	depth ++;

	/*
	 * The Borwein algo seems to always be faster than direct 
	 * summation, even when the direct-sum region is made quite 
	 * small, e.g. even when it is of radius less than 1/4.
	 * Never use direct summation.
	 */
#if 0
	if (0.0625>mod)
	{
		cpx_polylog_sum (plog, ess, zee, prec);
		return 0;
	}
#endif

	/* The zone of convergence for the Borwein algorithm is 
	 * |z^2/(z-1)| < 3.  If z is within this zone, then all is 
	 * well. If not, use the duplication formula to make 
	 * recursive calls, until the leaves of the recursion 
	 * are in this zone. 
	 *
	 * The algo seems to be more precise (!??) and have less 
	 * trouble when an even smaller bound is used, e.g.
	 * when |z^2/(z-1)| < 1.25. The non-recursive algo seems
	 * to choke up when it gets too close to z=1.
	 */

	/* den = | z^2/(z-1)|^2 */
	double den = polylog_get_zone (zre, zim);

	/* nterms = number of terms needed for computation */
	int nterms = polylog_terms_est (ess, zee, prec);

	/* To carry out the computation, an internal precision is needed 
	 * that is a bit higher than what the user asked for. This internal
	 * precision depends on the degree of the approximating polynomial.
	 * Thus, look at the available bits of precision, and decide if
	 * the calculation can be performed in that way. 
	 *
	 * The degree of internal precision available limits the largest
	 * effective order of the apprximating polynomial that can be used.
	 * Its pointless/erroneous to try to use a polynomial of degree
	 * more than "maxterms".
	 */
	int nbits = mpf_get_default_prec();
	int maxterms = nbits - (int) (3.321928095 *prec); /* log(10) / log(2) */

	// printf ("invoke-twrds, z=%g +i %g  den=%g nterms=%d, maxterms=%d\n", zre, zim, den, nterms, maxterms);

	/* If the z value is sufficently close to z=-1, then the Borwein 
	 * algorithm can be applied directly. So apply it. Oh, make sure
	 * that it doesn't take a hopeless numer of terms to get there.
	 */
	if ((den < 1.5) && (maxterms > nterms))
	{
		/* Use the larger, adjusted internal precision discussed above
		 * in the final calculation.
		 */
		prec += (int) (0.301029996 * nterms) +1;
		polylog_borwein (plog, ess, zee, nterms, prec);
		return 0;
	}

	/* Everything inside the unit circle is best handled by the
	 * recurse_away algorithm, which is proven, safe and effective.
	 * In fact, this should even be extended a bit, probably?
	 */
	if (mod <= 1.0)
	{
		rc = recurse_away_polylog (plog, ess, zee, prec, depth);
		return rc;
	}

	/* Use the polylog-hurwitz reflection formula, if the z value
	 * is sufficiently close to z=1. Basically, the hurwitz series 
	 * converges well when |q| < 0.5, or, in this case, if 
	 * |ln z / (2 pi i) | < 0.5. Note that flt point is good 
	 * enough for this test.
	 */
	if (log(mod) < 6.28)
	{
		rc = polylog_invert (plog, ess, zee, prec, depth);
		return rc;
	}

	/* If we are here, we are not close to either z=-1 or z=+1, and
	 * so the only option is to use the duplication formula to try 
	 * to get into one of these regions.
	 *
	 * Actually, with current algo, this statement should never be
	 * reached.
	 */
	printf ("splitsville-sqrt, z=%g +i %g  den=%g nterms=%d\n", zre, zim, den, nterms);

	rc = polylog_recurse_sqrt (plog, ess, zee, prec, depth);
	return rc;
}

int cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	int rc = recurse_towards_polylog (plog, ess, zee, prec, 0);
	if (rc)
	{
		cpx_set_ui (plog, 0,0);
		return rc;
	}
	return 0;
}

/* ============================================================= */
/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Caches intermediate results, so that overall performance is
 * considerably better if z is varied while s is held fixed.
 *
 * The magnitude of z must be less than one.
 */

void cpx_polylog_sum (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	int n;

	cpx_t s, z, zp, term;
	cpx_init (s);
	cpx_init (z);
	cpx_init (zp);
	cpx_init (term);

	cpx_set (s, ess);
	cpx_set (z, zee);
	cpx_set (zp, zee);
	cpx_set_ui (plog, 0, 0);

	/* Estimate the number of terms needed to sum over */
	double mag = polylog_modsq (zee);
	
	/* Domain error, should be less than one */
	if (1.0 <= mag)
	{
		fprintf (stderr, "cpx_polylog_sum(): Domain error, |z|=%g\n", sqrt(mag)); 
		return;
	}
	
	int nterms = -2.0 * prec *2.302585093 / log(mag);
	for (n=1; n<nterms; n++)
	{
		cpx_ui_pow_cache (term, n, s, prec);
		cpx_div (term, zp, term);

		cpx_add (plog, plog, term);
		cpx_mul (zp, zp, z);
	}

	cpx_clear (s);
	cpx_clear (z);
	cpx_clear (zp);
	cpx_clear (term);
}

/* ============================================================= */
/**
 * cpx_polylog_nint -- compute the polylogarithm at negetive integers
 *
 * At the negative integers, the polylog is a rational function,
 * meromorphic everywhere except for multiple poles at z=1.
 */

void cpx_polylog_nint (cpx_t plog, unsigned int negn, const cpx_t zee)
{
	int k;

	mpz_t stir, fac;
	mpz_init (stir);
	mpz_init (fac);

	cpx_t z, zp, term;
	cpx_init (z);
	cpx_init (zp);
	cpx_init (term);

	cpx_set (z, zee);
	cpx_sub_ui (zp, zee, 1, 0);
	cpx_recip (zp, zp);

	if (0 == negn)
	{
		cpx_mul (plog, z, zp);
		cpx_neg (plog, plog);
	}
	else
	{
		cpx_set_ui (plog, 0, 0);
		mpz_set_ui (fac, 1);
		cpx_set (z, zp);
		for (k=1; k<= negn+1; k++)
		{
			i_stirling_second (stir, negn+1, k);
			mpz_mul (stir, stir, fac);
			mpf_set_z (term[0].re, stir);
			mpf_set_ui (term[0].im, 0);

			cpx_mul (term, term, zp);
	
			cpx_add (plog, plog, term);
			cpx_mul (zp, zp, z);
			mpz_mul_ui (fac, fac, k);
		}

		if (0==negn%2)
		{
			cpx_neg (plog, plog);
		}
	}

	cpx_clear (z);
	cpx_clear (zp);
	cpx_clear (term);
	mpz_clear (stir);
	mpz_clear (fac);
}

/* ============================================================= */
/**
 * cpx_polylog_euler -- compute the polylogarithm from Hurwitz Euler.
 *
 * Combine two Hurwitz Euler-Maclaurin evaluations to obtain the polylogarithm.
 */

void cpx_polylog_euler (cpx_t zeta, const cpx_t ess, const cpx_t zee, int prec)
{
	int redo = 0;
	cpx_t q, tmp;
	cpx_init (q);
	cpx_init (tmp);

	/* Create a cache of commonly re-used values */
	static int cache_prec = -1;

	static mpf_t twopi, otp, log_twopi;
	static cpx_t phase, scale, cache_ess, s;

	if (-1 == cache_prec)
	{
		mpf_init (twopi);
		mpf_init (otp);
		mpf_init (log_twopi);

		cpx_init (phase);
		cpx_init (scale);
		cpx_init (s);
		cpx_init (cache_ess);
		cpx_set_ui(cache_ess, 123123123, 321321321);
	}

	if (cache_prec != prec)
	{
		redo = 1;
		cache_prec = prec;
		mpf_set_prec (twopi, 3.322*prec +50);
		mpf_set_prec (otp, 3.322*prec +50);
		mpf_set_prec (log_twopi, 3.322*prec +50);

		cpx_set_prec (phase, 3.322*prec +50);
		cpx_set_prec (scale, 3.322*prec +50);
		cpx_set_prec (s, 3.322*prec +50);
		cpx_set_prec (cache_ess, 3.322*prec +50);

		fp_two_pi (twopi, prec);

		/* otp = -1/2pi */
		mpf_set_ui (otp, 1);
		mpf_neg (otp, otp);
		mpf_div (otp, otp, twopi);

		fp_log (log_twopi, twopi, prec);
	}

	/* Recompute these values only if s differs from last time. */
	if(redo || !cpx_eq (ess, cache_ess, prec*3.322))
	{
		cpx_set (cache_ess, ess);
		cpx_ui_sub (s, 1, 0, ess);

		/* compute phase = e^{i pi s / 2} = (i)^s */
		cpx_times_mpf (tmp, s, twopi);
		cpx_div_ui (tmp, tmp, 4);
		cpx_times_i (tmp, tmp);
		cpx_exp (phase, tmp, prec);

		/* scale = gamma(s) / (2pi)^s */
		cpx_gamma_cache (scale, s, prec);
		cpx_times_mpf (tmp, s, log_twopi);
		cpx_neg (tmp, tmp);
		cpx_exp (tmp, tmp, prec);
		cpx_mul (scale, scale, tmp);
	}

	/* Compute q = ln z/(2pi i) */
	cpx_log (q, zee, prec);
	cpx_times_mpf (q, q, otp);
	cpx_times_i (q, q);

	/* exp (i pi s/2) * zeta (s,q) */
	cpx_hurwitz_euler (zeta, s, q, prec);
	cpx_mul (tmp, zeta, phase);

	/* exp (-i pi s/2) * zeta (s,1-q) */
	cpx_ui_sub (q, 1, 0, q);
	cpx_hurwitz_euler (zeta, s, q, prec);
	cpx_div (zeta, zeta, phase);
	cpx_add (zeta, zeta, tmp);

	cpx_mul (zeta, zeta, scale);

	cpx_clear (q);
	cpx_clear (tmp);
}

/* ============================================================= */
/**
 * cpx_periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 *
 * For 0<q<1/4 or for 3/4<q<1, this algo applies the duplication
 * formula, and calls itself recusrively, until 1/4<q<3/4, which
 * can then be evaluated at a single shot using teh Borwein algorithm.
 */
void cpx_periodic_zeta (cpx_t z, const cpx_t ess, const mpf_t que, int prec)
{
	mpf_t q, qf;
	mpf_init (q);
	mpf_init (qf);
	
	cpx_t s, sm;
	cpx_init (s);
	cpx_init (sm);

	mpf_set (q, que);
	mpf_floor (qf, q);
	mpf_sub (q, q, qf);

	cpx_set (s, ess);
	
	double fq = mpf_get_d (q);
	if ((1.0e-15 > fq) || (1.0e-15 > 1.0-fq))
	{
		// XXX should be more precise with the next order
		// q correction ... 
		// riemann_zeta (s.re, s.im, &z.re, &z.im);
		cpx_set_ui (z, 0, 0);
	}
	else if (mpf_cmp_d (q, 0.25) < 0) 
	{
		/* Use the duplication formula to get into convergent region */
		cpx_t ts, bt;
		cpx_init (ts);
		cpx_init (bt);

		/* sm = 1-s */
		cpx_neg (sm, s);
		cpx_add_ui (sm, sm, 1, 0);

		/* ts = 2^{1-s} */
		fp_log2 (qf, prec);
		cpx_times_mpf (sm, sm, qf);
		cpx_exp (ts, sm, prec);
		
		/* bt = pzeta (2q) * 2^{1-s} */
		mpf_mul_ui (qf, q, 2);
		cpx_periodic_zeta (bt, s, qf, prec);
		cpx_mul (bt, bt, ts);

		/* pzeta (q+0.5) */
		mpf_set_ui (qf, 1);
		mpf_div_ui (qf, qf, 2);
		mpf_add (qf, q, qf);
		cpx_periodic_zeta (z, s, qf, prec);
		cpx_sub (z, bt, z);
		
		cpx_clear (ts);
		cpx_clear (bt);
	}
	else if (mpf_cmp_d (q, 0.75) > 0) 
	{
		/* Use the duplication formula to get into convergent region */
		cpx_t ts, bt;
		cpx_init (ts);
		cpx_init (bt);

		/* sm = 1-s */
		cpx_neg (sm, s);
		cpx_add_ui (sm, sm, 1, 0);

		/* ts = 2^{1-s} */
		fp_log2 (qf, prec);
		cpx_times_mpf (sm, sm, qf);
		cpx_exp (ts, sm, prec);
		
		/* bt = pzeta (2q-1) * 2^{1-s} */
		mpf_mul_ui (qf, q, 2);
		mpf_sub_ui (qf, qf, 1);
		cpx_periodic_zeta (bt, s, qf, prec);
		cpx_mul (bt, bt, ts);

		/* pzeta (q-0.5) */
		mpf_set_ui (qf, 1);
		mpf_div_ui (qf, qf, 2);
		mpf_sub (qf, q, qf);
		cpx_periodic_zeta (z, s, qf, prec);
		cpx_sub (z, bt, z);

		cpx_clear (ts);
		cpx_clear (bt);
	}
	else
	{
		/* Normal case; within the convergence region */
		fp_two_pi (qf, prec);
		mpf_mul (qf, qf, q);

		fp_cosine (z[0].re, qf, prec);
		fp_sine (z[0].im, qf, prec);
		
		// cpx_polylog (z, s, z, prec);
		int nterms = polylog_terms_est (s, z, prec);
		if (4 < nterms)
		{
			polylog_borwein (z, s, z, nterms, prec);
		}
		else
		{
			fprintf (stderr, "Error: cpx_periodic_zeta() has bad terms estimate\n");
		}
	}
	
	mpf_clear (q);
	mpf_clear (qf);

	cpx_clear (s);
	cpx_clear (sm);
}

/* ============================================================= */
/**
 * cpx_periodic_beta -- Periodic beta function 
 *
 * Similar to periodic zeta, but with different normalization
 *
 * beta = 2 Gamma(s+1) (2\pi)^{-s} F(s,q)
 *
 * The implemented algorithm is to compute periodic zeta, and 
 * then renormalize and return the ressult.
 *
 * As of 22 December 2006, seems to be passing the tests -- 
 * that is, it gives the Bernoulli polynomials for integer s,
 * with all the right scale factors and signs, etc. Yay!
 */
void cpx_periodic_beta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	static cpx_t cache_s, scale;
	static int precision = 0;
	int redo = 0;

	if (!precision)
	{
		precision = prec;
		cpx_init (cache_s);
		cpx_init (scale);
	}

	if (precision < prec)
	{
		precision = prec;
		redo = 1;
		cpx_set_prec (cache_s, 3.322*prec+50);
		cpx_set_prec (scale, 3.322*prec+50);
	}

	if (redo || !cpx_eq (ess, cache_s, prec*3.322))
	{
		cpx_set (cache_s, ess);

		mpf_t two_pi;
		mpf_init (two_pi);

		cpx_t s, tps;
		cpx_init (s);
		cpx_init (tps);

		/* 2 gamma(s+1)/ (2pi)^s */
		cpx_add_ui (s, ess, 1,0);
		cpx_gamma_cache (scale, s, prec);

		/* times (2pi)^{-s} */
		fp_two_pi (two_pi, prec);
		cpx_neg (s, ess);
		cpx_mpf_pow (tps, two_pi, s, prec); 
		cpx_mul (scale, scale, tps);

		/* times two */
		cpx_times_ui (scale, scale, 2);
		cpx_clear (tps);
		cpx_clear (s);
		mpf_clear (two_pi);
	}

	cpx_periodic_zeta (zee, ess, que, prec);
	cpx_mul (zee, zee, scale);
}

/* ============================================================= */
/**
 * hurwitz_zeta -- Hurwitz zeta function
 *
 * Built up from the periodic zeta. Caches intermediate terms, and so
 * performance is much better if s is held const, while q is varied.
 * Expects the input value of q to be between 0 and 1.
 *
 * "built up" means periodic zeta is computed (twice) and then summed
 * with appropriate factors, to compute hurwitz zeta.
 */

static void hurwitz_zeta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	static cpx_t cache_s, piss, niss, scale;
	static int precision = 0;
	int redo = 0;

	if (!precision)
	{
		precision = prec;
		cpx_init (cache_s);
		cpx_init (piss);
		cpx_init (niss);
		cpx_init (scale);
	}
	if (precision < prec)
	{
		precision = prec;
		redo = 1;
		cpx_set_prec (cache_s, 3.322*prec+50);
		cpx_set_prec (piss, 3.322*prec+50);
		cpx_set_prec (niss, 3.322*prec+50);
		cpx_set_prec (scale, 3.322*prec+50);
	}

	mpf_t t;
	mpf_init (t);

	cpx_t s, zm;
	cpx_init (s);
	cpx_init (zm);

	/* s = 1-ess */
	cpx_neg (s, ess);
	cpx_add_ui (s, s, 1, 0);

	if (redo || !cpx_eq (s, cache_s, prec*3.322))
	{
		cpx_set (cache_s, s);

		cpx_t tps;
		cpx_init (tps);

		/* exp (i pi s/2) */
		fp_pi_half (t, prec);
		cpx_times_mpf (piss, s, t);
		cpx_times_i (piss, piss);
		cpx_exp (piss, piss, prec);
		cpx_recip (niss, piss);

		/* gamma(s)/ (2pi)^s */
		cpx_gamma_cache (scale, s, prec);

		/* times (2pi)^{-s} */
		fp_two_pi (t, prec);
		cpx_neg (s, s);
		cpx_mpf_pow (tps, t, s, prec); 
		cpx_neg (s, s);
		cpx_mul (scale, scale, tps);

		/* times two */
		cpx_clear (tps);
	}
	
	/* F(s,q) and F(s, 1-q) */
	cpx_periodic_zeta (zee, s, que, prec);
	mpf_ui_sub (t, 1, que);
	cpx_periodic_zeta (zm, s, t, prec);

	/* assemble the thing */
	cpx_mul (zm, zm, piss);
	cpx_mul (zee, zee, niss);
	cpx_add (zee, zee, zm);
	cpx_mul (zee, zee, scale);

	cpx_clear (s);
	cpx_clear (zm);
	mpf_clear (t);
}

/* ============================================================= */
/**
 * cpx_hurwitz_zeta -- Hurwitz zeta function
 *
 * Built up from the periodic zeta. Caches intermediate terms, and so
 * performance is much better if s is held const, while q is varied.
 * The value for q must be positive; the algo gets slow if q is large.
 */

void cpx_hurwitz_zeta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	cpx_t s, term;
	mpf_t q;
	cpx_init (s);
	cpx_init (term);
	mpf_init (q);
	cpx_set (s, ess);
	mpf_set (q, que);

	/* Make sure q is between 0 and 1, by subtracting the integer part */
	long nq = mpf_get_si (q);
	mpf_sub_ui (q, q, nq);

	hurwitz_zeta (zee, s, q, prec);

	cpx_neg (s, s);
	int k;
	for (k=0; k<nq; k++)
	{
		/* subtract off the leading terms -- the 1/(q+k)^s */
		cpx_mpf_pow (term, q, s, prec);
		cpx_sub (zee, zee, term);
		mpf_add_ui (q, q, 1);
	}

	cpx_clear (s);
	cpx_clear (term);
	mpf_clear (q);
}

/* ============================================================= */
/**
 * cpx_hurwitz_taylor -- Hurwitz zeta function Taylor series
 *
 * Implement the Hurwitz zeta as a taylor expansion about q=0
 * (pulling out the leading 1/q^s term to handle uniquely)
 */

void cpx_hurwitz_taylor (cpx_t zee, const cpx_t ess, const cpx_t que, int prec)
{
	cpx_t s, sn, q, qn, bin, term;
	cpx_init (s);
	cpx_init (sn);
	cpx_init (q);
	cpx_init (qn);
	cpx_init (bin);
	cpx_init (term);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm, aterm;
	mpf_init (maxterm);
	mpf_init (aterm);
	fp_epsilon (maxterm, 2*prec);

	cpx_set (s, ess);
	cpx_set (q, que);

	/* Compute 1/q^s if the real part of q is near 0. 
	 * Else navigate the waters of the Hurwitz branch cut.
	 */
	cpx_neg (s, s);
	cpx_set_ui (zee, 0, 0);
	while (mpf_cmp_d (q[0].re, 0.5) < 0)
	{
		cpx_pow (qn, q, s, prec);
		cpx_add (zee, zee, qn);
		mpf_add_ui (q[0].re, q[0].re, 1);
	}
	mpf_sub_ui (q[0].re, q[0].re, 1);
	while (mpf_cmp_d (q[0].re, 1.5) > 0)
	{
		mpf_sub_ui (q[0].re, q[0].re, 1);
		cpx_pow (qn, q, s, prec);
		cpx_sub (zee, zee, qn);
	}
	if (mpf_cmp_d (q[0].re, 0.5) > 0)
	{
		mpf_sub_ui (q[0].re, q[0].re, 1);
	}
	cpx_neg (s, s);

double qre = mpf_get_d (q[0].re);
double qim = mpf_get_d (q[0].im);
double mod = sqrt(qre*qre+qim*qim);
if (0.9 < mod) {
printf ("oh no mr bill! q=%g +i%g m=%g\n", qre, qim, mod);
goto punt;
}
	/* Now do a power series in -q */
	cpx_neg (q, q);
	cpx_set_ui (qn, 1, 0);
	cpx_sub_ui (sn, s, 1, 0);
	int n = 0;
	while (1)
	{
		/* s+n-1 */
		/* The caching version uses precomputed values,
		 * making te algo faster. This also means that 
		 * sn does not need to be incremented. */
		// cpx_binomial (bin, sn, n);
		// cpx_add_ui (sn, sn, 1, 0);
		cpx_binomial_sum_cache (bin, sn, n);

		/* zeta_cache returns vale at s+n */
		cpx_borwein_zeta_cache (term, s, n, prec);
		// cpx_borwein_zeta (term, sn, prec);
		cpx_mul (term, term, bin);
		cpx_mul (term, term, qn);
		cpx_add (zee, zee, term);

		/* Don't go no farther than this */
		cpx_mod_sq (aterm, term);
		if (mpf_cmp (aterm, maxterm) < 0) break;

		n++;

		/* qn = q^n */
		cpx_mul (qn, qn, q);
	}

punt:
	cpx_clear (s);
	cpx_clear (sn);
	cpx_clear (q);
	cpx_clear (qn);
	cpx_clear (bin);
	cpx_clear (term);
	mpf_clear (aterm);
	mpf_clear (maxterm);
}

/* =========================================================== */
/**
 * cpx_hurwitz_euler -- Hurwitz zeta function via Euler-Maclaurin algo
 *
 * This function computes the value of the Hurwitz zeta function
 * using an Euler-Maclaurin summation to obtain an estimate.
 *
 * The algorithm appears to work.
 */

static void zeta_euler_fp(cpx_t zeta, cpx_t ess, mpf_t q, int em, int prec)
{
	int k;
	cpx_t s, spoch, term, deriv;
	cpx_init (s);
	cpx_init (spoch);
	cpx_init (term);
	cpx_init (deriv);

	cpx_neg (s, ess);

	cpx_set_ui (zeta, 0, 0);
	/* sum over 1/(k+q)^s  from k=0 to k=M-1 */
	for (k=0; k<em; k++)
	{
		fp_pow_rc (term, k, q, s, prec);
		cpx_add (zeta, zeta, term);
	}

	/* deriv = 1/(M+q)^s */
	fp_pow_rc (deriv, em, q, s, prec);
	
	/* Add another (1/2) of 1 /(M+q)^s */
	cpx_div_ui (term, deriv, 2);
	cpx_add (zeta, zeta, term);

	mpf_t eps, fact, emq, ft;
	mpf_init (eps);
	mpf_init (fact);
	mpf_init (emq);
	mpf_init (ft);

	mpq_t bern;
	mpq_init (bern);
	
	/* emq = M+q */
	mpf_add_ui (emq, q, em);

	/* term = 1/(s-1)*(q+M)^{s-1} */
	cpx_times_mpf (term, deriv, emq);
	cpx_add_ui (s, s, 1, 0);
	cpx_div (term, term, s);
	cpx_sub (zeta, zeta, term);

	/* emq = 1/(M+q) */
	mpf_ui_div (emq, 1, emq);
	cpx_times_mpf (deriv, deriv, emq);

	/* emq = 1/(M+q)^2 */
	mpf_mul (emq, emq, emq);
	
	mpf_set_ui (fact, 1);
	mpf_div_ui (fact, fact, 2);
	
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_set (spoch, s);

	fp_epsilon (eps, 2*prec);
	
	k = 1;
	while (1)
	{
		/* ft = B_2k / (2k)! */
		q_bernoulli (bern, 2*k);
		mpf_set_q (ft, bern);
		mpf_mul (ft, ft, fact);
		
		cpx_times_mpf (term, deriv, ft);
		cpx_mul (term, term, spoch);
		cpx_add (zeta, zeta, term);

		cpx_mod_sq (ft, term);
		if (mpf_cmp (ft, eps) < 0) break;
#if 0
		double t = mpf_get_d (ft);
		printf ("M=%d Q=%d bern=%g ", em, k, t);
#endif

		k++;
		
		mpf_div_ui (fact, fact, (2*k-1)*2*k);
		cpx_times_mpf (deriv, deriv, emq);
		mpf_add_ui (s[0].re, s[0].re, 1);
		cpx_mul (spoch, spoch, s);
		mpf_add_ui (s[0].re, s[0].re, 1);
		cpx_mul (spoch, spoch, s);
	}
	
	mpq_clear (bern);
	mpf_clear (fact);
	mpf_clear (emq);
	mpf_clear (ft);
	cpx_clear (s);
	cpx_clear (spoch);
	cpx_clear (term);
	cpx_clear (deriv);
} 

static void zeta_euler(cpx_t zeta, cpx_t ess, cpx_t q, int em, int prec)
{
	int k;
	cpx_t s, emq, spoch, term, deriv;
	cpx_init (s);
	cpx_init (emq);
	cpx_init (spoch);
	cpx_init (term);
	cpx_init (deriv);

	cpx_neg (s, ess);
	cpx_set (emq, q);

	cpx_set_ui (zeta, 0, 0);
	/* sum over 1/(k+q)^s  from k=0 to k=M-1 */
	for (k=0; k<em; k++)
	{
		cpx_pow_rc (term, k, emq, s, prec);
		cpx_add (zeta, zeta, term);
	}

	/* deriv = 1/(M+q)^s */
	cpx_pow_rc (deriv, em, emq, s, prec);
	
	/* Add another (1/2) of 1 /(M+q)^s */
	cpx_div_ui (term, deriv, 2);
	cpx_add (zeta, zeta, term);

	mpf_t eps, fact, ft;
	mpf_init (eps);
	mpf_init (fact);
	mpf_init (ft);

	mpq_t bern;
	mpq_init (bern);
	
	/* emq = M+q */
	mpf_add_ui (emq[0].re, emq[0].re, em);

	/* term = 1/(s-1)*(q+M)^{s-1} */
	cpx_mul (term, deriv, emq);
	cpx_add_ui (s, s, 1, 0);
	cpx_div (term, term, s);
	cpx_sub (zeta, zeta, term);

	/* emq = 1/(M+q) */
	cpx_recip (emq, emq);
	cpx_mul (deriv, deriv, emq);

	/* emq = 1/(M+q)^2 */
	cpx_mul (emq, emq, emq);
	
	mpf_set_ui (fact, 1);
	mpf_div_ui (fact, fact, 2);
	
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_set (spoch, s);

	fp_epsilon (eps, 2*prec);
	
	k = 1;
	while (1)
	{
		/* ft = B_2k / (2k)! */
		q_bernoulli (bern, 2*k);
		mpf_set_q (ft, bern);
		mpf_mul (ft, ft, fact);
		
		cpx_times_mpf (term, deriv, ft);
		cpx_mul (term, term, spoch);
		cpx_add (zeta, zeta, term);

		cpx_mod_sq (ft, term);
		if (mpf_cmp (ft, eps) < 0) break;
#if 0
		double t = mpf_get_d (ft);
		printf ("M=%d Q=%d bern=%g ", em, k, t);
#endif

		k++;
		
		mpf_div_ui (fact, fact, (2*k-1)*2*k);
		cpx_mul (deriv, deriv, emq);
		mpf_add_ui (s[0].re, s[0].re, 1);
		cpx_mul (spoch, spoch, s);
		mpf_add_ui (s[0].re, s[0].re, 1);
		cpx_mul (spoch, spoch, s);
	}
	
	mpq_clear (bern);
	mpf_clear (fact);
	mpf_clear (ft);
	cpx_clear (s);
	cpx_clear (spoch);
	cpx_clear (term);
	cpx_clear (deriv);
} 

void cpx_hurwitz_euler_fp(cpx_t zeta, cpx_t ess, mpf_t q, int prec)
{
	/* really really really bad estimates to the bounds */
	// int em = prec + 12;
	int em = prec/2 + 5;

	zeta_euler_fp (zeta, ess, q, em, prec);
}

void cpx_hurwitz_euler(cpx_t zeta, cpx_t ess, cpx_t q, int prec)
{
	/* really really really bad estimates to the bounds */
	int em = prec + 12;

	zeta_euler (zeta, ess, q, em, prec);
}

/* ============================================================= */
/* Implement algorithm 1 from Cohen, Villegas, Zagier et al 
 * The naive implementation below is a total failure, I don't know why.
 * Clearly, for real s, one might expect failure, since the series
 * is not alternating. However, for s with a large imaginary component,
 * the terms in the series are oscillatory, i.e. quasi-alternating,
 * and, for whatever reason, doesn't work there.
 */

#ifdef BORKEN_DOESNT_WORK_DONT_KNOW_WHY
void cpx_pade_hurwitz_zeta (cpx_t hur, const cpx_t ess, const mpf_t que, int prec)
{
	int k;
	int nterms;

	nterms = 1.31*prec;
	
	cpx_t s, term;
	cpx_init (s);
	cpx_init (term);
	cpx_set (s, ess);

	mpf_t b,c,d,q;
	mpf_init (q);
	mpf_init (b);
	mpf_init (c);
	mpf_init (d);
	mpf_set (q, que);

	/* d = (3+sqrt(8))^n */
	mpf_set_ui (d, 8);
	mpf_sqrt (d, d);
	mpf_add_ui (d,d,3);
	mpf_pow_ui (d,d,nterms);

	/* d = 0.5 * (d + 1/d) */
	mpf_ui_div (c, 1, d);
	mpf_add (d, d, c);
	mpf_div_ui (d, d, 2);

	mpf_neg (c, d);
	mpf_set_ui (b, 1);
	mpf_neg (b, b);

	cpx_set_ui (hur, 0,0);
	for (k=0; k<nterms-1; k++)
	{
		mpf_sub (c, b, c);
		fp_pow_rc (term, k, q, s, prec);
		cpx_times_mpf (term, term, c);
		cpx_add (hur, hur, term);

		mpf_mul_ui (b,b, 2*(k+nterms)*(k-nterms));
		mpf_div_ui (b,b, (2*k+1)*(k+1));
	}

	cpx_div_mpf (hur, hur, d);
	
	mpf_clear (b);
	mpf_clear (c);
	mpf_clear (d);
	mpf_clear (q);
	cpx_clear (s);
	cpx_clear (term);
}
#endif

/* ============================================================= */

/** 
 * test_bernoulli_poly - compare periodic zeta to the Bernoulli poly's
 *
 * The Bernoulli polynomials have a fourier transform that is the 
 * periodic zeta function. 
 *
 * Test is now passing with flying colors
 */
#if 0
int test_bernoulli_poly (int n)
{
	cplex zl, zh;

	cplex s, z;
	s.im = 0.0;
	s.re = n;
	double q;
	for (q=-0.2; q<=1.2; q+=0.02)
	{
		// zl = periodic_zeta (s, q);
		// zh = periodic_zeta (s, 1.0-q);
		zl = periodic_beta (s, q);
		zh = periodic_beta (s, 1.0-q);
		if (n%2) {
			z = cplex_sub (zl,zh);
		} else {
			z = cplex_add (zl,zh);
		}
		
		double bs;
		if (0 == n%2)
		{
	 		bs = z.re;
			if (n%4 == 0) bs = -bs;
		}
		if (n%2)
		{
			bs = -z.im;
			if (n%4 == 3) bs = -bs;
		}

		// bs *= factorial (n) * pow (2.0*M_PI, -n);
		bs *= 0.5;

		/* short-circuit the above, test directly */
		cplex ess;
		ess.re = 1-n;
		ess.im = 0;
		cplex hz = hurwitz_zeta(ess,q);
		bs = -n * hz.re;
		
		// double b = q*q-q+1.0/6.0;
		double b = bernoulli_poly (n,q);
		
		printf ("q=%5.3g	bs=%13.10g	bernoulli=%13.10g	reldiff=%6.3g\n", q, bs, b, (bs-b)/b);
	}
}

/* ============================================================= */
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

main (int argc, char * argv[])
{
	int n;
	double en;

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <n>\n", argv[0]);
		exit (1);
	}
	n = atoi (argv[1]);
	en = atof (argv[1]);

	// test_zeta_values (n);
	// test_bernoulli_poly (n);

// #define HURWITZ_ZETA
#ifdef HURWITZ_ZETA
	/* As of 22 December 2006, this test seems to be passing, 
	 * with decent accuracy, for anything with real part less than about 8
	 */
	cplex s;
	s.im = 0.0;
	double q=0.5;
	for (s.re = 1.05; s.re < n; s.re += 0.1)
	{
		cplex hz= hurwitz_zeta (s, q);
		
		double zeta = gsl_sf_hzeta (s.re, q);
		
		printf ("s=%5.3g	algo=%12.10g	exact=%12.10g	reldiff=%6.3g\n", s.re, hz.re, zeta, (hz.re-zeta)/zeta);
	}
#endif

}
#endif

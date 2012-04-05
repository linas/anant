/* polylog-bug.c -- Ad-hoc bug report -- Test calc of Dilogarithm. */

/* Last change: 27-Mar-2012.  Ken Roberts  */

/* Calc Li_2(u+iv)
 * then Li_2(u+iv)
 * again.  Expect consistent results.
 * Looking for possible bug re caching
 * of info within cpx_polylog logic.
 *
 * This exhibits a caching bug in anant-0.2.0   
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <db_185.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gmp.h>
#include "mp-complex.h"
#include "mp-polylog.h"

ulong gprec;
ulong aprec;

int main (
  int    argc,
  char * argv[])
{

  double u, v;
  cpx_t s, z, pl1b, pl1e;
  int  res;

  gprec = 1000;
  mpf_set_default_prec (gprec);
  aprec = 50;

  cpx_init (s);
  cpx_init (z);

  cpx_set_ui (s, 2, 0);

  u = 0.000005;
  v = 1.377128;

  cpx_set_d (z, u, v);

  printf ("\n");
  gmp_printf ("s = %10Ff + i %10Ff\n",
              s[0].re, s[0].im);
  gmp_printf ("z = %10Ff + i %10Ff\n",
              z[0].re, z[0].im);

  cpx_init (pl1b);
  cpx_init (pl1e);

  res = cpx_polylog (pl1b, s, z, aprec);
  gmp_printf ("Rtn %d, Li_2(z) = %10Ff + i %10Ff\n",
              res, pl1b[0].re, pl1b[0].im);

  res = cpx_polylog (pl1e, s, z, aprec);
  gmp_printf ("Rtn %d, Li_2(z) = %10Ff + i %10Ff\n",
              res, pl1e[0].re, pl1e[0].im);

  printf ("\n");
  exit (0);
};

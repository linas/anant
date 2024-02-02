
Anant  -- Algorithmic 'n Analytic Number Theory
===============================================

<p style="text-align: center;"><b>
<i>Version 0.3.0</i><br>
Linas Vepstas February 2024<br>
linasvepstas@gmail.com
</b></p>


This project contains ad-hoc implementations of assorted analytic
functions of interest in number theory, including the gamma function,
the Riemann zeta function, the polylogarithm, and the Minkowski
question mark function. The implementation uses the Gnu Multi-Precision
library (GMP) to perform all low-level operations.  The code herein
is licensed under the terms of the Gnu GPLv3 license.

This project is *NOT* meant to be a replacement for other, more
established multi-precision systems, such as MPFR or PARI/GP.  It is
meant to be a staging area for implementations of functions that have
not (yet!) received much attention in the more established packages.
Users are strongly encouraged to port the contents of this package
to other systems. Please! 
Seriously: the stuff here is still cutting-edge, fifteen years on.
I really would like to see the "serious" math packages actually
implement these functions!

A note about floating-point exceptions: many of the special functions
computed here have poles at various values.  These will show up as
mysterious floating-point exceptions deep within the code.  If you
get an exception, make sure you are not evaluating a function at a pole!

This package has its origins as a collection of tools & utilities for
the benefit of the author.  As such, it was never really intended for
public consumption, and thus, will not have the usual amenities of
established projects, such as clear documentation, a website, unit
test cases, or even a robust build system.  Caveat Emptor! It works
great for me, and I use it all the time; but I am me, and not you.

There are several publications that describe this code, or make use of
it. The most notable is this:

Vepstas, L. (2008) ["An efficient algorithm for accelerating
the convergence of oscillatory series, useful for computing the
polylogarithm and Hurwitz zeta functions"]
(http://link.springer.com/article/10.1007/s11075-007-9153-8).
*Numerical Algorithms* vol. **47** issue 3: pp. 211-252.
[arXiv:math.CA/0702243](https://arxiv.org/abs/math/0702243).
doi:10.1007/s11075-007-9153-8.

The most fun thing created with this code is this:
[Linas' Art Gallery](https://linas.org/art-gallery).

The most verbose thing created with this code is this:
[The Modular Group and Fractals](https://linas.org/math/sl2z.html)
<i>An exploration of the relationship between Fractals, the Riemann Zeta,
the Modular Group Gamma, the Farey Fractions and the Minkowski Question
Mark.</i>

Supported functions
-------------------
Many of the functions below are not "difficult"; what makes the code
here unique is that many of these use caching and partial computation
to avoid repeated computations. In some cases, this caching allows the
functions to be particularly fast when called in "sequential" order,
as would naturally occur in summations. This can provide speedups of
over a thousand-fold, when used correctly!

Another, apparently "unique" feature, as compared to established
libraries, is that most of the functions are handled on the complex
plane.  Thus, the system comes with basic complex-number support.

Some, but not all of the code is protected by spinlocks, and so can be
used in parallel, i.e. in a multi-threaded setting. This can also give
a big boost in performance.  This code is used regularly on a 24-core
CPU. It works.

Arbitrary precision constants
-----------------------------
* sqrt(3)/2, log(2)
* e, e^pi
* pi, 2pi, pi/2, sqrt(2pi), log(2pi), 2/pi
* Euler-Mascheroni const
* Riemann zeta(1/2)

Combinatorial functions:
------------------------
* Rising pochhammer symbol (integer)
* Partition function (integer)
* Reciprocal factorial
* Sequential binomial coefficient
* Stirling Numbers of the First Kind
* Stirling Numbers of the Second Kind
* Bernoulli Numbers
* Binomial transform of power sum
* Rising pochhammer symbol (real)    i.e. (s)_n for real s
* Rising pochhammer symbol (complex) i.e. (s)_n for complex s
* Binomial coefficient (complex)     i.e. (s choose n) for complex s
* High-speed, cached sequential-access binomial coefficient.

Elementary functions:
---------------------
* pow, exp, log, sine, cosine, tangent for real, complex arguments
* arctan, arctan2 for real argument
* log(1-x) for real, complex x
* sqrt for complex argument

Classical functions:
--------------------
* gamma (factorial) for real, complex argument
* polylogarithm, using multiple algorithms: Borwein-style, Euler-Maclaurin
* polylogarithm on multiple sheets (monodromy)
* Periodic zeta function
* Hurwitz zeta function, using multiple algorithms; complex arguments.
* Riemann zeta function, using multiple algorithms: Borwein, Hasse, brute-force
  for integer, real and complex arguments.
* Confluent hypergeometric function, complex arguments

Number-theoretic functions:
---------------------------
* General complex-valued harmonic number
* Gauss-Kuzmin-Wirsing operator matrix elements
* Minkowski Question Mark function (Stern-Brocot tree), and its inverse
* Taylor's series coefficients for the topologist's sin -- sin(2pi/(1+x))

Utilities:
----------
* Fill in values of completely multiplicative arithmetic function.
* Ordinary and exponential generating functions for arithmetic functions.
* Euler re-summation of convergent series (speeds convergence).
* Newton series interpolation of arithmetic series.
* Powell's method for zero-finding on complex plane (noise-cancelling variant).
* Root isolation on complex plane using Sagraloff-Yap (2011) algorithm.


Pre-requisites, Compiling, Installing, Testing
----------------------------------------------
This package requires a copy of the Berkeley DB database to be
available. The database is used to cache certain intermediate values, to
improve performance of various internal algorithms.

This package has minimal build support. `cd` to the src directory, and
`make`.  If you want to install the files somewhere, you will have to do
this by hand, or custom-tailor to suit your needs.

There is a unit test, rather ad-hoc in nature, and it is not
"user-friendly".  It will report some errors in the last few decimal
places of various routines, depending on how it was invoked. It is up
to you to figure out if these are serious errors or not.  Caveat Emptor!

(I mean, its 'error-free' for me, i.e. 'good enough'.  There is nothing
in here that is horribly broken, as far as I know.  If quibbling over
the last few decimal places is important to you, you might have a
different opinion. If you're reasonably careful, and actually think
about what you are actually doing, things will go well.)

Patches to improve the build system (and anything else that annoys you)
are gladly accepted.

Precision
---------
Most of the algorithms deal with precision issues in a fairly ad-hoc
kind of way.  Many/most routines require an argument specifying the
number of decimal places of desired precision, and will typically
return answers that are accurate from about 90% to 100% of the specified
precision. However, many of the algorithms use internal, intermediate
results that need to be maintained at a higher level of precision than
the "desired" precision. Thus, correct usage requires that the user
specify an mpf_set_default_prec() that is 10% to 50% larger than the
desired precision of the results.  The proper amount to use is up to
you to figure out!  A reasonable rule-of-thumb seems to be to use
mpf_set_default_prec(5*desired_decimal_places) -- noting that log_2(10)
is 2.3.

Again -- this may strike you as hacky; its good enough for what I need.
Patches to improve the state of things are accepted.

Example Usage
-------------
The below provides an example of how to use the functions in this
library.

      // Standard include headers
      #include <gmp.h>
      #include <stdio.h>
      #include "mp-polylog.h"
      #include "mp-misc.h"

      int main()
      {
         cpx_t plog, ess, zee;  // Complex variant of mpf_t
         int nbits;
         int decimal_prec;

         // decimal_prec is the number of decimal places of desired
         // precision.
         //
         // 3.3 is equal to log_2(10), and is used to convert decimal
         // places to number of binary bits.
         //
         // The +600 adds some extra "padding precision" for
         // intermediate calculations. Most algorithms require some
         // fair amount of additional bits of precision to be used in
         // computing intermediate results.  The precise amount needed
         // is somewhat ad-hoc, and not well-characterized for the
         // different functions; typically, an extra 20% to 50% is
         // needed.
         //
         decimal_prec = 500;
         nbits = 3.3*decimal_prec + 600;
         mpf_set_default_prec (nbits);

         // Initialization
         cpx_init (plog);
         cpx_init (ess);
         cpx_init (zee);

         // Set values for which computation will be done.
         cpx_set_d(ess, 2.1, 0.0);
         cpx_set_d(zee, 0.5, 0.0);

         // Compute ...
         int rc = cpx_polylog(plog, ess, zee, decimal_prec);

         // Check for error conditions
         if (0 != rc)
         {
             printf("Error occurred during computation! rc=%d\n", rc);
             return 1;
         }
         cpx_prt("Answer is ", plog);
         printf("\n");

         return 0;
      }


Current repository:
-------------------
in Git, on Github:
   https://github.com/linas/anant

Older versions (2005 through 2012) can be found in Bazaar, on Launchpad:
   https://launchpad.net/anant

Source tarballs are available there too.

Related projects:
-----------------
* [C library for arbitrary-precision interval arithmetic](http://arblib.org/)
  from Fredrik Johansson.  This library includes Airy and Bessel
  functions, the error function, as well as the modular j-function and
  others. [Github repo](https://github.com/fredrik-johansson/arb/)

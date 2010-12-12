/**
 * mp-zerofind.h
 * Locate zeros of a function. 
 *
 * Linas Vepstas December 2010
 */

#include "mp-complex.h"


/* =============================================== */
/**
 * find_zero.
 * Numerically locate the zero of a complex-valued function.
 * 
 * @func function whose zeros are to be found.
 *       func takes z as input, returns f as output. 
 *       'nprec' is the suggested decimal precisiton at which 'fun'
 *       should perform its calculations.
 * @initial_z initial suggestion for the location of the zero.
 * @e1, @e2 initial suggestions for a bounding ellipse. These are
 *       taken to be two vectors, specifying the major and minor axes of
 *       an ellipse, centered at 'initial_z'. The true zero is presumed
 *       to lie inside of, or at least, close to, this initial ellipse.
 * @ndigits number of decimal digits of accuracy to which the zero
 *       should be searched for.
 * @nprec number of digits of decimal precision to which intermediate
 *       terms will be maintained.
 *
 * @returns 0 if result is valid, else an error code.
 *
 * This implements Powell's method, slightly adapted; the adaptations
 * are meant to imporve convergence when 'func' is extremely noisy, 
 * i.e. when any sort of quardatic behaviour is obscured by
 * high-freqency noise.
 *
 * A description of Powell's method can be found in Press, Teukolsky,
 * Vetterling, Flannery, "Numerical Recipes in C, 2nd ed.", Cambridge 
 * U Press, 1999.
 */

int find_zero(cpx_t result,
              void (*func)(cpx_t f, cpx_t z, int nprec),
              cpx_t initial_z,
              cpx_t e1, cpx_t e2,
              int ndigits, int nprec);

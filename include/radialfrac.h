//---------------------------------------------------------------------------
//
// $Id: radialfrac.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------

/* radialfrac.h
 * evaluate the continued fraction */

#ifndef RADIALFRAC_H_SEEN
#define RADIALFRAC_H_SEEN

#include "teukolskydefs.h"

/* estimate of the precision of the results of the fractions, as
 * a multiple of their term sizes
 * This is generally a very generous bound, because we want to catch outliers
 */
#define TERM_NOISE 1e-14

/* evaluate the fraction with increasing n
 * - a_0 / (b1 - a_1 / (... */
Real     plusfrac(Real    nu, Real epsilon, Real q, int m, Real lambda,
		  Real *term);
Complex cplusfrac(Complex nu, Real epsilon, Real q, int m, Real lambda,
		  Real *term);

/* evaluate the fraction with decreasing n
 * - a_-1 / (b_-1 - a_-2 / (... */
Real     minusfrac(Real    nu, Real epsilon, Real q, int m, Real lambda,
		   Real *term);
Complex cminusfrac(Complex nu, Real epsilon, Real q, int m, Real lambda,
		   Real *term);

/* evaluate the whole thing */
Real     radialfrac(Real    nu, Real epsilon, Real q, int m, Real lambda,
		    Real *term);
Complex cradialfrac(Complex nu, Real epsilon, Real q, int m, Real lambda,
		    Real *term);

/* Evaluate the continued fraction on the special line Re(nu) = -1/2
 * It is guaranteed to be real there */
Real radialfrac_half(Real imnu, Real epsilon, Real q, int m, Real lambda,
		     Real *term);

#endif /* RADIALFRAC_H_SEEN */

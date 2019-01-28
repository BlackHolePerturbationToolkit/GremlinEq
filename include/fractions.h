//---------------------------------------------------------------------------
//
// $Id: fractions.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
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
/* renangmom.h
 * compute the renormalized angular momentum nu by solving the radial
 * continued fraction equation
 */

#ifndef RENANGMOM_H_SEEN
#define RENANGMOM_H_SEEN

#include "teukolskydefs.h"

/* calculates nu
 * returns 1 if it suceeds, 0 else */
int renangmom(Real epsilon, Real q, int l, int m, Real lambda, Complex *nu);

/* calulates nu assuming it is real.  returns 1 if it suceeds, 0 else */
int renangmom_real(Real epsilon, Real q, int l, int m, Real lambda, Real *nu);

/* uses about the closest thing to a brute force root finder that I could
 * come up with to search for real nu. Returns 1 if it succeeds */
int renangmom_real_divide_search(Real epsilon, Real q, int l, int m,
				 Real lambda, Real *nu);

/* calulates nu assuming it is on the special line.
 * Returns 1 if it suceeds, 0 else */
int renangmom_half(Real epsilon, Real q, int l, int m, Real lambda,
		   Real *imnu);

/* calculates nu assuming it is on the line Re(nu) = 1 (so any int is
 * then a solution)
 * Returns 1 if it suceeds, 0 else */
int renangmom_iint(Real epsilon, Real q, int l, int m, Real lambda,
		   Real *imnu);

/* Calculates nu when epsilon is too large for the region checks to work
 * Assumes nu is complex */
int renangmom_far(Real epsilon, Real q, int l, int m, Real lambda,
		  Complex *nu);

/* a very rough guess at Im(nu), used internally.
 * returns the imaginary part of the nu such that beta(0) vanishes */
Real nu_predictor(Real epsilon,Real q,int m,Real lambda);

#endif /* RENANGMOM_H_SEEN */
/* specialradial.h
 * Functions to calculate properties of the radial continued fraction
 * at special points
 */

#ifndef SPECIALRADIAL_H_SEEN
#define SPECIALRADIAL_H_SEEN

#include "teukolskydefs.h"

/* estimate of the precision of the results of the special functions, as
 * a multiple of their term sizes */
#define SPECIAL_TERM_NOISE 3e-15

/* calculate the value of the function at nu = -1/2
 * If term is not NULL, the address it points to will be set to the magnitude
 * of one of the terms in the final addition */
Real radial_half_value(Real epsilon, Real q, int m, Real lambda, Real *term);

/* calculate the slope of the function at nu = 1
 * If term is not NULL, the address it points to will be set to the magnitude
 * of one of the terms in the final addition */
Real radial_int_slope(Real epsilon, Real q, int m, Real lambda, Real *term);

#endif /* SPECIALRADIAL_H_SEEN */

//---------------------------------------------------------------------------
//
// $Id: renangmom.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------

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

//---------------------------------------------------------------------------
//
// $Id: fsum.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------

/* fsum.h
 * Perform sums of the expansion coefficients of the radial Teukolsky equation
 */

#ifndef FSUM_H_SEEN
#define FSUM_H_SEEN

#include "teukolskydefs.h"

/* Calculate the asymptitic amplitudes
 * Any of the pointers can be null */
void asympt_amps(Complex nu, Real epsilon, Real q, int m, Real lambda,
		 Complex *b_trans, Complex *b_inc, Complex *b_ref,
		 Complex *c_trans);

/* perform a straight sum of the f's from -\infty to \infty (to some
 * approximation), normalizing f_0 = 1 */
Complex fsum(Complex nu, Real epsilon, Real q, int m, Real lambda);

/* Compute the K_\nu factor */
Complex kfactor(Complex nu, Real epsilon, Real q, int m, Real lambda);

/* compute A_{-}^\nu */
Complex aminus(Complex nu, Real epsilon, Real q, int m, Real lambda);

/* compute R_0^\nu
 * if deriv is non-NULL it will be set to d/dr R_0 */
Complex rzero(Complex nu, Real epsilon, Real q, int m, Real lambda, Real x,
	      Complex *deriv);

/* compute R^{in}_{lmw} from hypergeometric functions
 * if deriv is non-NULL it will be set to d/dr R^{in} */
Complex rin_hyper(Complex nu, Real epsilon, Real q, int m, Real lambda,
		  Real x, Complex *deriv);

/* compute R_C */
Complex rcoulomb(Complex nu, Real epsilon, Real q, int m, Real lambda, Real z);

/* compute R^{in}_{lmw} from Coulomb wave functions */
Complex rin_coulomb(Complex nu, Real epsilon, Real q, int m, Real lambda,
		    Real z);

/* compute R^{in}_{lmw}
 * if deriv is non-NULL it will be set to d/dr R^{in} */
Complex rin(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
	    Complex *deriv);

/* compute R^{up}_{lmw} using Tricomi functions
 * if deriv is non-NULL it will be set to d/dr R^{up} */
Complex rup_tricomi(Complex nu, Real epsilon, Real q, int m, Real lambda,
		    Real z, Complex *deriv);

/* compute R^{up}_{lmw}
 * if deriv is non-NULL it will be set to d/dr R^{up} */
Complex rup(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
	    Complex *deriv);

/* compute R^{in}_{lmw} from hypergeometric functions optimized for
 * the case where r is small
 * if deriv is non-NULL it will be set to d/dr R^{in} */
Complex rin_small(Complex nu, Real epsilon, Real q, int m, Real lambda,
		  Real x, Complex *deriv);

#endif /* !FSUM_H_SEEN */

//---------------------------------------------------------------------------
//
// $Id: specialradial.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------


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

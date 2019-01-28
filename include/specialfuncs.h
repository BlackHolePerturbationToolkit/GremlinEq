//---------------------------------------------------------------------------
//
// $Id: specialfuncs.h,v 1.1 2013/05/17 15:30:12 sahughes Exp $
//
//---------------------------------------------------------------------------

/* funcs.h
 * headers for the special functions routines
 */

#ifndef FUNCS_H_SEEN
#define FUNCS_H_SEEN

#include <complex>

/* These functions have been declared using the built in types because some
 * of them (mainly gammln) will not give enough precision for a long double */

/* logarithm of the gamma function */
std::complex<double> gammln(const std::complex<double> x);

/* logarithm of the sine function */
std::complex<double> sinln(const std::complex<double> x);

/* hypergeometric function 2F1(n1,n2;d1;x)
 * only works for abs(x) < 1
 * it has some issues for large magnitude n1,n2,d1, particularly near x = -1
 */
std::complex<double>
hypergeom2F1(std::complex<double> n1, std::complex<double> n2,
	     std::complex<double> d1, std::complex<double> x);

/* confluent hypergeometric function 1F1(n1;d1;x) */
std::complex<double>
hypergeom1F1(std::complex<double> n1, std::complex<double> d1,
	     std::complex<double> x);

/* irregular Tricomi hypergeometric function U(a;b;x) */
std::complex<double>
hypergeomU(std::complex<double> a, std::complex<double> b,
	   std::complex<double> x);

/* irregular Tricomi hypergeometric function U(a;b;x)*Gamma(a-b+1)/Gamma(1-b)
 */
std::complex<double>
hypergeomU_reduced(std::complex<double> a, std::complex<double> b,
		   std::complex<double> x);

#endif /* !FUNCS_H_SEEN */

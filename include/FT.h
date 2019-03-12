//---------------------------------------------------------------------------
//
// $Id: FT.h,v 1.2 2018/08/02 20:19:01 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// A class for interacting with the Teukolsky equation solver described
// by Fujita and Tagoshi
//

#ifndef FT_H_SEEN
#define FT_H_SEEN

#include "teukolskydefs.h"

#define FT_MAXDIVS 10
#define TEST_LENGTH 4

#define FT_PRFX FT::

//! Fujita Tagoshi Class
/*! FT defines methods for computing solutions to the homogeneous Teukolsky equation as laid out in papers by Ryuichi Fujita and Hideyushi Tagoshi */
class FT {
public:
  FT(const int l, const int m, const Real r, const Real a, const Real omega,
     const Real lambda, const Real tolerance);
  FT(const int l, const int m, const Real p, const Real ecc, const Real a,
     const Real omega, const Real lambda, const Real tolerance);
  
  Complex Bin()    { return b_inc; } /*!< asymptotic amplitude */
  Complex Btrans() { return b_trans; } /*!< asymptotic amplitudes */
  Complex Ctrans() { return c_trans; } /*!< asymptotic amplitudes */

  void CalcRFields(const Real rad, const int DoH); /*!< set the point at which the fields are calculated */

  Complex TeukRin()       { return rteuk_in; } /*!< value of the homogeneous solutions */
  Complex TeukRup()     { return rteuk_up; } /*!< value of the homogeneous solutions */
  Complex dr_TeukRin()    { return drteuk_in; } /*!< value of the homogeneous solutions */
  Complex dr_TeukRup()  { return drteuk_up; } /*!< value of the homogeneous solutions */
  Complex ddr_TeukRin()   { return ddrteuk_in; } /*!< value of the homogeneous solutions */
  Complex ddr_TeukRup() { return ddrteuk_up; } /*!< value of the homogeneous solutions */

  Real Accuracy_in()   { return accuracy_in; } /*!< the actual error bound */
  Real Accuracy_up() { return accuracy_up; } /*!< the actual error bound */

  Real request_precision_in  (const Real p) {
    Real tmp = rin_request_precision;
    rin_request_precision = Fmax(p,REAL_EPSILON);
    return tmp;
  } /*!< requesting precision */

  Real request_precision_up(const Real p) {
    Real tmp = rup_request_precision;
    rup_request_precision = Fmax(p,REAL_EPSILON);
    return tmp;
  } /*!< requesting precision */
  Real request_precision_in  () { return rin_request_precision; } /*!< requesting precision */
  Real request_precision_up() { return rup_request_precision; } /*!< requesting precision */

  int terms_evaluated() {
    int tmp = hypergeom_terms;
    hypergeom_terms = 0;
    return tmp;
  }

  void set_gmp(int flag);
  int get_gmp() { return gmp_on; }

  
  inline Real K(const Real rad) {
    return((rad*rad + a*a)*omega - m*a);
  } /*!< some quantities in the Teukolsky equation */

  inline Real dr_K(const Real rad) {
    return(2.*rad*omega);
  } /*!< some quantities in the Teukolsky equation */

  inline Real ddr_K() {
    return(2.*omega);
  } /*!< some quantities in the Teukolsky equation */

  inline Complex teuk_potential(const Real rad) {
    return - K(rad) / (rad*rad - 2*rad + a*a) * (K(rad) + 4*(rad-1)*I)
      + 8 * omega * rad * I + lambda;
  } /*!< some quantities in the Teukolsky equation */

  inline Complex get_ddrteuk(const Real rad,
			     const Complex rteuk, const Complex drteuk); /*!< some quantities in the Teukolsky equation */

  int choose_solvers(const Real low, const Real high,
		     const Real *lowerrs, const Real *higherrs,
		     const int isin, const int oldchoice,
		     Real **locptr, int **choiceptr, const int *maxchoiceptr); /*!< some quantities in the Teukolsky equation */

  void geterrs(const Real r, Real *errs, int isin);

  Complex callsolver(const int isin, const int which, const Real r,
		     const Real epsilon, const Complex nu,
		     const Real precision, Complex *deriv);

  Complex callin(const Real r, Complex *deriv);

  Complex callup(const Real r, Complex *deriv);

  void getlocs(const int isin);

  int l, m;
  Real a, p, e, omega, lambda;
  Real tolerance, accuracy_in, accuracy_up;
  Real rin_request_precision, rup_request_precision;

  Complex nu;
  Complex b_trans, b_inc, c_trans;

  Complex rteuk_in, rteuk_up;
  Complex drteuk_in, drteuk_up;
  Complex ddrteuk_in, ddrteuk_up;

  int up_choices[FT_MAXDIVS];
  Real up_choice_locs[FT_MAXDIVS];
  int in_choices[FT_MAXDIVS];
  Real in_choice_locs[FT_MAXDIVS];

  Complex test_renangmoms[TEST_LENGTH];

  int gmp_on;
  int hypergeom_terms;

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

/*!< evaluate the fraction with increasing n
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

/* calculates the value of the fraction at -1/2 and the slope at 1.
 * Returns 1 if the values appear to have the right sign */
int get_deciders(Real epsilon, Real q, int m, Real lambda,
		 Real *halfval, Real *oneslope);

/* calculates nu assuming it is real by calling the other
 * renangmom_real* functions.  returns 1 on success, 0 else */
int renangmom_real_search_procedure(Real epsilon, Real q, int l,
				    int m, Real lambda, Complex *nu);

/* assuming the root is near an integer or half integer, fit the two
 * approaches to epsilon with polynomials and extrapolate.
 * returns 1 on success, 0 else */
int singularity_fit(Real epsilon, Real q, int l, int m, Real lambda,
		    int guessint, int guessint2, Real spacing, Complex *nu);

/* calulates nu assuming it is real.  returns 1 if it suceeds, 0 else */
int renangmom_real(Real epsilon, Real q, int l, int m, Real lambda,
		   Real detect_limit, Real *nu);

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

/* search for the root assuming the function is well behaved near guess and
 * that the root is closer to guess than to a half integer */
int renangmom_real_guess(Real epsilon, Real q, int m, Real lambda, Real guess,
			 Real *nu);

/* approximate the root using a polynomial approximation in large nu */
Real fractions_poly_approx(Real epsilon, Real q, int m, Real lambda);

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
#define SPECIAL_TERM_NOISE 1e-12

/* calculate the value of the function at nu = -1/2
 * If term is not NULL, the address it points to will be set to the magnitude
 * of one of the terms in the final addition */
Real radial_half_value(Real epsilon, Real q, int m, Real lambda, Real *term);

/* calculate the slope of the function at nu = 1
 * If term is not NULL, the address it points to will be set to the magnitude
 * of one of the terms in the final addition */
Real radial_int_slope(Real epsilon, Real q, int m, Real lambda, Real *term);

#endif /* SPECIALRADIAL_H_SEEN */

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
		  Real x, Real precision, Complex *deriv);

/* compute R_C */
Complex rcoulomb(Complex nu, Real epsilon, Real q, int m, Real lambda, Real z,
		 Complex *deriv);

/* compute R^{in}_{lmw} from Coulomb wave functions */
Complex rin_coulomb(Complex nu, Real epsilon, Real q, int m, Real lambda,
		    Real z, Real precision, Complex *deriv);

/* compute R^{in}_{lmw}
 * if deriv is non-NULL it will be set to d/dr R^{in} */
/*Complex rin(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
  Complex *deriv);*/

/* compute R^{up}_{lmw} using Tricomi functions
 * if deriv is non-NULL it will be set to d/dr R^{up} */
Complex rup_tricomi(Complex nu, Real epsilon, Real q, int m, Real lambda,
		    Real z, Real precision, Complex *deriv);

/* compute R^{up}_{lmw}
 * if deriv is non-NULL it will be set to d/dr R^{up} */
/*Complex rup(Complex nu, Real epsilon, Real q, int m, Real lambda, Real r,
  Complex *deriv);*/

/* compute R^{in}_{lmw} from hypergeometric functions optimized for
 * the case where r is small
 * if deriv is non-NULL it will be set to d/dr R^{in} */
Complex rin_small(Complex nu, Real epsilon, Real q, int m, Real lambda,
		  Real x, Real precision, Complex *deriv);

/* compute R^{up}_{lmw} using Hypergeometric functions
 * if deriv is non-NULL it will be set to d/dr R^{up} */
Complex rup_hyper(Complex nu, Real epsilon, Real q, int m, Real lambda,
		  Real x, Real precision, Complex *deriv);

#endif /* !FSUM_H_SEEN */

/* funcs.h
 * headers for the special functions routines
 */

#ifndef FUNCS_H_SEEN
#define FUNCS_H_SEEN

#include <complex>

/* USE_GMP
 * undef             = gmp disabled at build time
 * defined but false = gmp forced
 * defined and true  = gmp controlable at runtime
 */
#ifndef NO_GMP
#  define USE_GMP 1
#endif

#define HYPERGEOM_REAL_TYPE long double

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
std::complex<HYPERGEOM_REAL_TYPE>
hypergeom2F1(std::complex<double> n1, std::complex<double> n2,
	     std::complex<double> d1, std::complex<double> x);

/* confluent hypergeometric function 1F1(n1;d1;x) */
std::complex<HYPERGEOM_REAL_TYPE>
hypergeom1F1(std::complex<double> n1, std::complex<double> d1,
	     std::complex<double> x);

/* irregular Tricomi hypergeometric function U(a;b;x) */
std::complex<HYPERGEOM_REAL_TYPE>
hypergeomU(std::complex<double> a, std::complex<double> b,
	   std::complex<double> x);

/* irregular Tricomi hypergeometric function U(a;b;x)*Gamma(a-b+1)/Gamma(1-b)
 */
std::complex<HYPERGEOM_REAL_TYPE>
hypergeomU_reduced(std::complex<double> a, std::complex<double> b,
		   std::complex<double> x);

#ifdef USE_GMP
std::complex<HYPERGEOM_REAL_TYPE>
hypergeom2F1_gmp(std::complex<double> n1, std::complex<double> n2,
		 std::complex<double> d1, std::complex<double> x,
		 int precision);
std::complex<HYPERGEOM_REAL_TYPE>
hypergeom1F1_gmp(std::complex<double> n1, std::complex<double> d1,
		 std::complex<double> x, int precision);
#endif /* USE_GMP */

#endif /* !FUNCS_H_SEEN */

};

#endif /* !FT_H_SEEN */

//---------------------------------------------------------------------------
//
// $Id: hypergeom.cc,v 1.5 2014/04/08 19:36:44 sahughes Exp $
//
//---------------------------------------------------------------------------

/* hypergeom.cc
 * functions for computing hypergeometric functions
 */

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "funcs.h"

/* toggling gmp uses the object interface */
#  ifdef USE_GMP
#    undef USE_GMP
#    define USE_GMP 0
#  endif
#endif

#include <cstdio>
#include <cfloat>
#include "complex_ops.h"

#define hpmax(a,b)				\
  ({						\
    typeof(a) _a = (a);				\
    typeof(b) _b = (b);				\
    _a > _b ? _a : _b;				\
  })

#ifdef USE_GMP
#include <gmpxx.h>

/* minimum estimated precision in hyper* before the gmp version is called */
#define MIN_PRECISION 1e-14

/* how much precision over the estimated need to use */
#define GMP_EXTRA_PRECISION 10 /* 40 */

#endif /* USE_GMP */

std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeom2F1(std::complex<double> n1, std::complex<double> n2,
		     std::complex<double> d1, std::complex<double> x) {
  if(abs(x) >= 1 && x.imag() != 0.) {
    fprintf(stderr,"Argument %g%+gi too large in hypergeom2F1\n",
	    real(x),imag(x));
    throw("Argument too large in hypergeom2F1");
  }

  int n;
  std::complex<HYPERGEOM_REAL_TYPE> sum = 1;
  std::complex<HYPERGEOM_REAL_TYPE> term = 1;

#ifdef USE_GMP
  HYPERGEOM_REAL_TYPE maxabs = 1;
  HYPERGEOM_REAL_TYPE abssum;
#endif

  for(n=0;;n++) {
    term *= (n1 + n)*(n2 + n) / ((d1 + n) * (n + 1)) * x;
    if(sum == sum + term) {
#ifdef HYPERGEOM_TRACK_TERMS
      hypergeom_terms += n+1;
#endif /* HYPERGEOM_TRACK_TERMS */
#ifdef USE_GMP
#  if USE_GMP
      if(gmp_on) {
#  endif
	abssum = abs(sum);
	if(abssum * MIN_PRECISION < maxabs * DBL_EPSILON)
	  return hypergeom2F1_gmp(n1, n2, d1, x,
				  (int) log2(maxabs / abssum)
				  + DBL_MANT_DIG + GMP_EXTRA_PRECISION);
#  if USE_GMP
      }
#  endif
#endif
      return sum;
    }
    sum += term;
    if(isnan(real(sum)) || isnan(imag(sum)))
      return std::complex<HYPERGEOM_REAL_TYPE>(NAN,NAN);
#ifdef USE_GMP
#  if USE_GMP
    if(gmp_on) {
#  endif
      abssum = abs(sum);
      if(abssum > maxabs) maxabs = abssum;
#  if USE_GMP
    }
#  endif
#endif
  }
}

std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeom1F1(std::complex<double> n1, std::complex<double> d1,
		     std::complex<double> x) {
  int n;
  std::complex<HYPERGEOM_REAL_TYPE> sum = 1;
  std::complex<HYPERGEOM_REAL_TYPE> term = 1;

#ifdef USE_GMP
  HYPERGEOM_REAL_TYPE maxabs = 1;
  HYPERGEOM_REAL_TYPE abssum;
#endif

  for(n=0;;n++) {
    term *= (n1 + n) / ((d1 + n) * (n + 1)) * x;
    if(sum == sum + term) {
#ifdef HYPERGEOM_TRACK_TERMS
      hypergeom_terms += n+1;
#endif /* HYPERGEOM_TRACK_TERMS */
#ifdef USE_GMP
#  if USE_GMP
      if(gmp_on) {
#  endif
	abssum = abs(sum);
	if(abssum * MIN_PRECISION < maxabs * DBL_EPSILON)
	  return hypergeom1F1_gmp(n1, d1, x,
				  (int) log2(maxabs / abssum)
				  + DBL_MANT_DIG + GMP_EXTRA_PRECISION);
#  if USE_GMP
      }
#  endif
#endif
      return sum;
    }
    sum += term;
    if(isnan(real(sum)) || isnan(imag(sum)))
      return std::complex<HYPERGEOM_REAL_TYPE>(NAN,NAN);
#ifdef USE_GMP
#  if USE_GMP
    if(gmp_on) {
#  endif
      abssum = abs(sum);
      if(abssum > maxabs) maxabs = abssum;
#  if USE_GMP
    }
#  endif
#endif
  }
}

/* it may be better to reimplement this using the method described in
 * Allasia and Besenghi (1987), but this is simpler */
std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeomU(std::complex<double> a, std::complex<double> b,
		   std::complex<double> x) {
  return
    static_cast<std::complex<HYPERGEOM_REAL_TYPE> >
    (exp(gammln(b-1.0) - gammln(a) + (1.0-b)*log(x)))
    * hypergeom1F1(a-b+1.0,2.0-b,x)
    + static_cast<std::complex<HYPERGEOM_REAL_TYPE> >
    (exp(gammln(1.0-b) - gammln(a-b+1.0)))
    * hypergeom1F1(a,b,x);
}

std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeomU_reduced(std::complex<double> a, std::complex<double> b,
			   std::complex<double> x) {
  return
    static_cast<std::complex<HYPERGEOM_REAL_TYPE> >
    (exp(gammln(b-1.0) + gammln(a-b+1.0) - gammln(1.0-b) - gammln(a)
	 + (1.0-b)*log(x)))
    * hypergeom1F1(a-b+1.0,2.0-b,x)
    + hypergeom1F1(a,b,x);
}

#ifdef USE_GMP
/* I have no idea why the result of these calculations are meaningful.
 * I would expect the error from the large terms in the sum to be very
 * high, but something conspires to make the errors cancel and give
 * a relatively small derivative.
 * (In 2F1 at -13.7-250.4I, -11.7-2I, -27.4, 0.15 all derivatives are
 * well under 100 ulps/ulp, according to Mathematica)
 */
std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeom2F1_gmp(std::complex<double> n1, std::complex<double> n2,
			 std::complex<double> d1, std::complex<double> x,
			 int precision) {
  const mpf_class n1_r(real(n1),precision);
  const mpf_class n1_i(imag(n1),precision);
  const mpf_class n2_r(real(n2),precision);
  const mpf_class n2_i(imag(n2),precision);
  const mpf_class d1_r(real(d1),precision);
  const mpf_class d1_i(imag(d1),precision);
  const mpf_class x_r(real(x),precision);
  const mpf_class x_i(imag(x),precision);

  mpf_class sum_r(1,precision);
  mpf_class sum_i(0,precision);
  mpf_class term_r(1,precision);
  mpf_class term_i(0,precision);

  mpf_class mult_r(0,precision);
  mpf_class mult_i(0,precision);

  mpf_class tmp(0,precision);

  int n;
  HYPERGEOM_REAL_TYPE maxabs = 1;
  HYPERGEOM_REAL_TYPE abssum;

  for(n=0;;n++) {
    mult_r = x_r*(n1_r + n) - x_i*n1_i;
    mult_i = x_r*n1_i + x_i*(n1_r + n);

    tmp = mult_r*(n2_r + n) - mult_i*n2_i;
    mult_i = mult_r*n2_i + mult_i*(n2_r + n);
    mult_r = tmp;

    tmp = mult_r*(d1_r + n) + mult_i*d1_i;
    mult_i = - mult_r*d1_i + mult_i*(d1_r + n);
    mult_r = tmp;

    tmp = d1_r + n;
    tmp = tmp*tmp;
    tmp += d1_i*d1_i;
    tmp *= n + 1;

    mult_r /= tmp;
    mult_i /= tmp;

    tmp = term_r*mult_r - term_i*mult_i;
    term_i = term_r*mult_i + term_i*mult_r;
    term_r = tmp;

    if(hpmax(abs(sum_r),abs(sum_i)) * DBL_EPSILON >
       hpmax(abs(term_r),abs(term_i))) {
#ifdef HYPERGEOM_TRACK_TERMS
      hypergeom_terms += n+1;
#endif /* HYPERGEOM_TRACK_TERMS */
      abssum = hypot(sum_r.get_d(),sum_i.get_d());
      if(abssum * exp2(precision) * MIN_PRECISION < maxabs)
	return hypergeom2F1_gmp(n1, n2, d1, x,
				(int) log2(maxabs / abssum)
				+ DBL_MANT_DIG + GMP_EXTRA_PRECISION);
      return std::complex<HYPERGEOM_REAL_TYPE>(sum_r.get_d(), sum_i.get_d());
    }
    sum_r += term_r;
    sum_i += term_i;
    abssum = hypot(sum_r.get_d(),sum_i.get_d());
    if(abssum > maxabs) maxabs = abssum;
  }
}

std::complex<HYPERGEOM_REAL_TYPE>
FT_PRFX hypergeom1F1_gmp(std::complex<double> n1, std::complex<double> d1,
			 std::complex<double> x, int precision) {
  const mpf_class n1_r(real(n1),precision);
  const mpf_class n1_i(imag(n1),precision);
  const mpf_class d1_r(real(d1),precision);
  const mpf_class d1_i(imag(d1),precision);
  const mpf_class x_r(real(x),precision);
  const mpf_class x_i(imag(x),precision);

  mpf_class sum_r(1,precision);
  mpf_class sum_i(0,precision);
  mpf_class term_r(1,precision);
  mpf_class term_i(0,precision);

  mpf_class mult_r(0,precision);
  mpf_class mult_i(0,precision);

  mpf_class tmp(0,precision);

  int n;
  HYPERGEOM_REAL_TYPE maxabs = 1;
  HYPERGEOM_REAL_TYPE abssum;

  for(n=0;;n++) {
    mult_r = x_r*(n1_r + n) - x_i*n1_i;
    mult_i = x_r*n1_i + x_i*(n1_r + n);

    tmp = mult_r*(d1_r + n) + mult_i*d1_i;
    mult_i = - mult_r*d1_i + mult_i*(d1_r + n);
    mult_r = tmp;

    tmp = d1_r + n;
    tmp = tmp*tmp;
    tmp += d1_i*d1_i;
    tmp *= n + 1;

    mult_r /= tmp;
    mult_i /= tmp;

    tmp = term_r*mult_r - term_i*mult_i;
    term_i = term_r*mult_i + term_i*mult_r;
    term_r = tmp;

    if(hpmax(abs(sum_r),abs(sum_i)) * DBL_EPSILON >
       hpmax(abs(term_r),abs(term_i))) {
#ifdef HYPERGEOM_TRACK_TERMS
      hypergeom_terms += n+1;
#endif /* HYPERGEOM_TRACK_TERMS */
      abssum = hypot(sum_r.get_d(),sum_i.get_d());
      if(abssum * exp2(precision) * MIN_PRECISION < maxabs)
	return hypergeom1F1_gmp(n1, d1, x,
				(int) log2(maxabs / abssum)
				+ DBL_MANT_DIG + GMP_EXTRA_PRECISION);
      return std::complex<HYPERGEOM_REAL_TYPE>(sum_r.get_d(), sum_i.get_d());
    }
    sum_r += term_r;
    sum_i += term_i;
    abssum = hypot(sum_r.get_d(),sum_i.get_d());
    if(abssum > maxabs) maxabs = abssum;
  }
}
#endif /* USE_GMP */

//---------------------------------------------------------------------------
//
// $Id: SWSHSpherical.cc,v 1.9 2019/01/26 14:38:10 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// The functions in this file are used to compute spin weighted
// spherical harmonics.  They are used by the class SWSH.
//
// Scott Hughes, 24 July 1998
//

#include <cmath>
#include "Globals.h"
#include "SWSH.h"

//
// The prefactor ("A(l, m)" in my notes) that appears before the
// spin-weight zero spherical harmonics.
//

Real SWSH::pref_A_func(const int l, const int m)
{
  Real ans;

  ans = sqrt((2.*((Real)l) + 1.)/(4.*M_PI));
  ans *= exp(0.5*(gsl_sf_lngamma(l - m + 1) - gsl_sf_lngamma(l + m + 1)));
  return ans;
}

//
// Recurrence relation for associated Legendre polynomial of
// x = cos(theta).
//

Real SWSH::plgndr(const int l, const int m, const Real x)
{
  Real fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if (m < 0 || m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine plgndr" << endl;
    exit(0);
  }

  pmm = 1.; // The value for m = 0.
  if (m > 0) {
    somx2 = sqrt((1. - x)*(1. + x));
    fact = 1.;
    for (i = 1; i <= m; i++) {
      pmm *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1 = x*(2*m + 1)*pmm;
    if (l == (m + 1))
      return pmmp1;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	pll = (x*(2*ll - 1)*pmmp1 - (ll + m - 1)*pmm)/(ll - m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}

//
// Recurrence relation for associated Legendre polynomial of
// x = cos(theta) divided by sqrt(1 - x^2) = sin(theta).  Only
// allows l >= 1 since this routine is used for the spherical
// harmonics of spin weight abs(s) > 0.
//

Real SWSH::plgndr_over_sin(const int l, const int m, const Real x)
{
  Real fact, pll_over_sin, pmm_over_sin, pmmp1_over_sin, somx2;
  int i, ll;

  if (l < 1 || m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine plgndr_over_sin" << endl;
    exit(0);
  }

  if (m == 0)
    pmm_over_sin = 1./sqrt((1. - x)*(1. + x));
  else pmm_over_sin = -1;  // The value for m = 1.
  if (m > 1) {
    somx2 = sqrt((1. - x)*(1. + x));
    fact = 3.;
    for (i = 2; i <= m; i++) {
      pmm_over_sin *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return pmm_over_sin;
  else {
    pmmp1_over_sin = x*(2*m + 1)*pmm_over_sin;
    if (l == (m + 1))
      return pmmp1_over_sin;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	pll_over_sin = (x*(2*ll - 1)*pmmp1_over_sin - 
		      (ll + m - 1)*pmm_over_sin)/(ll - m);
	pmm_over_sin = pmmp1_over_sin;
	pmmp1_over_sin = pll_over_sin;
      }
      return pll_over_sin;
    }
  }
}

//
// Recurrence relation for associated Legendre polynomial of
// x = cos(theta) divided by 1 - x^2 = sin^2(theta).  Only
// allows l >= 2 since this routine is used for the spherical
// harmonics of spin weight abs(s) > 1.
//

Real SWSH::plgndr_over_sinsqr(const int l, const int m, const Real x)
{
  Real fact, pll_over_sinsqr, pmm_over_sinsqr, pmmp1_over_sinsqr,
    somx2;
  int i, ll;

  if (l < 2 || m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine plgndr_over_sinsqr" << endl;
    exit(0);
  }

  if (m == 0) pmm_over_sinsqr = 1./((1. - x)*(1. + x));
  else if (m == 1) pmm_over_sinsqr = -1./sqrt((1. - x)*(1. + x));
  else pmm_over_sinsqr = 3.; // The value for m = 2.
  if (m > 2) {
    somx2 = sqrt((1. - x)*(1. + x));
    fact = 5.;
    for (i = 3; i <= m; i++) {
      pmm_over_sinsqr *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return pmm_over_sinsqr;
  else {
    pmmp1_over_sinsqr = x*(2*m + 1)*pmm_over_sinsqr;
    if (l == (m + 1))
      return pmmp1_over_sinsqr;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	pll_over_sinsqr = (x*(2*ll - 1)*pmmp1_over_sinsqr - 
		      (ll + m - 1)*pmm_over_sinsqr)/(ll - m);
	pmm_over_sinsqr = pmmp1_over_sinsqr;
	pmmp1_over_sinsqr = pll_over_sinsqr;
      }
      return pll_over_sinsqr;
    }
  }
}

//
// Recurrence relation for the derivative with respect to
// x = cos(theta) of the associated Legendre polynomial
// times sqrt(1 - x^2) = sin(theta).  Only allows l >= 1 since this
// routine is used for the spherical harmonics of spin weight
// abs(s) > 0.
//

Real SWSH::sin_times_dplgndr(const int l, const int m, const Real x)
{
  Real fact, sin_times_dpll, sin_times_dpmm,
    sin_times_dpmmp1, somx2;
  int i, ll;

  if (l < 1 || m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine sin_times_dplgndr" << endl;
    exit(0);
  }

  sin_times_dpmm = m*x; // The value for m = 1, times m (great for m = 0...)
  somx2 = sqrt((1. - x)*(1. + x));
  if (m > 1) {
    fact = 3.;
    for (i = 2; i <= m; i++) {
      sin_times_dpmm *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return sin_times_dpmm;
  else {
    sin_times_dpmmp1 = (2*m + 1)*(somx2*plgndr(m, m, x) + x*sin_times_dpmm);
    if (l == (m + 1))
      return sin_times_dpmmp1;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	sin_times_dpll = ((2*ll - 1)*(x*sin_times_dpmmp1 +
				    somx2*plgndr(ll - 1, m, x)) -
			  (ll + m - 1)*sin_times_dpmm)/(ll - m);
	sin_times_dpmm = sin_times_dpmmp1;
	sin_times_dpmmp1 = sin_times_dpll;
      }
      return sin_times_dpll;
    }
  }
}

//
// Recurrence relation for the derivative with respect to
// x = cos(theta) of the associated Legendre polynomial.
//

Real SWSH::dplgndr(const int l, const int m, const Real x)
{
  Real fact, dpll, dpmm, dpmmp1, somx2;
  int i, ll;

  if (m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine dplgndr" << endl;
    exit(0);
  }

  if (l == 0) return(0.);

  somx2 = sqrt((1. - x)*(1. + x));
  if (m == 1) dpmm = x/somx2;
  else dpmm = -3.*m*x; // The value for m = 2 (times m/2; great for m = 0...)
  if (m > 2) {
    fact = 5.;
    for (i = 3; i <= m; i++) {
      dpmm *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return dpmm;
  else {
    dpmmp1 = (2*m + 1)*(plgndr(m, m, x) + x*dpmm);
    if (l == (m + 1))
      return dpmmp1;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	dpll = ((2*ll - 1)*(x*dpmmp1 + plgndr(ll - 1, m, x)) -
			  (ll + m - 1)*dpmm)/(ll - m);
	dpmm = dpmmp1;
	dpmmp1 = dpll;
      }
      return dpll;
    }
  }
}

//
// Recurrence relation for the second derivative with respect
// to x = cos(theta) of the associated Legendre polynomial times
// (1 - x^2) = sin^2(theta).  Only allows l >= 2 since this routine
// is used for the spherical harmonics of spin weight abs(s) > 0.
//

Real SWSH::sinsqr_times_ddplgndr(const int l, const int m, const Real x)
{
  Real fact, sinsqr_times_ddpll, sinsqr_times_ddpmm,
    sinsqr_times_ddpmmp1, somx2;
  int i, ll;

  if (l < 2 || m > l || fabs(x) > 1.) {
    cerr << "Bad arguments in routine sinsqr_times_ddplgndr" << endl;
    exit(0);
  }

  somx2 = sqrt((1. - x)*(1. + x));
  if (m == 1)
    sinsqr_times_ddpmm = 1./somx2;
  else
    sinsqr_times_ddpmm = m*3.*(x*x*(m - 1) - 1.);
  // The value for m = 2, modulo factors (great for m = 0...)
  if (m > 2) {
    fact = 5.;
    for (i = 3; i <= m; i++) {
      sinsqr_times_ddpmm *= -fact*somx2;
      fact += 2.;
    }
  }
  if (l == m)
    return sinsqr_times_ddpmm;
  else {
    sinsqr_times_ddpmmp1 = (2*m + 1)*(2.*somx2*somx2*dplgndr(m, m, x) +
				    x*sinsqr_times_ddpmm);
    if (l == (m + 1))
      return sinsqr_times_ddpmmp1;
    else {
      for (ll = m + 2; ll <= l; ll++) {
	sinsqr_times_ddpll = ((2*ll - 1)*(x*sinsqr_times_ddpmmp1 +
					2.*somx2*somx2*dplgndr(ll - 1, m, x)) -
			      (ll + m - 1)*sinsqr_times_ddpmm)/(ll - m);
	sinsqr_times_ddpmm = sinsqr_times_ddpmmp1;
	sinsqr_times_ddpmmp1 = sinsqr_times_ddpll;
      }
      return sinsqr_times_ddpll;
    }
  }
}

//
// The below five functions calculate the three varieties of spherical
// harmonic that are needed.  Note that the recurrence relations for the
// various Legendre polynomials above are valid only for m non-negative;
// to handle negative m, we simply need to call the Legendre routines
// with abs(m) instead of m.  The various sums of those Legendre functions
// that have functional dependence on m are done with the "real" (not
// absolute valued) m.
//
// Note I am leaving off the factor exp(i m phi), since that phi dependence
// is separated away in the Teukolsky function.
//
Real SWSH::pos2Y(const int l, const int m, const Real x)
{
  Real pref, sinsqr_times_ddlegendreP, dlegendreP,
    legendreP_over_sinsqr;
  Real dl; // l recast as a Real - avoid integer overflows
  Real sum;

  dl = (Real)l;
  pref = pref_A_func(l, abs(m))/sqrt((dl - 1.)*dl*(dl + 1.)*(dl + 2.));
  if (m < 0 && m%2) pref *= -1;
  sinsqr_times_ddlegendreP = sinsqr_times_ddplgndr(l, abs(m), x);
  if (abs(m) > 0) {
    dlegendreP = dplgndr(l, abs(m), x);
    legendreP_over_sinsqr = plgndr_over_sinsqr(l, abs(m), x);
    sum = sinsqr_times_ddlegendreP + 2.*m*dlegendreP +
      (((Real)(m*m)) + 2.*m*x)*legendreP_over_sinsqr;
  } else sum = sinsqr_times_ddlegendreP;

  return pref*sum;
}

Real SWSH::pos1Y(const int l, const int m, const Real x)
{
  Real pref, sin_times_dlegendreP, legendreP_over_sin;
  Real dl; // l recast as a Real - avoid integer overflows
  Real sum;

  dl = (Real)l;
  pref = pref_A_func(l, abs(m))/sqrt(dl*(dl + 1.));
  if (m < 0 && m%2) pref *= -1;
  sin_times_dlegendreP = sin_times_dplgndr(l, abs(m), x);
  if (abs(m) > 0) {
    legendreP_over_sin = plgndr_over_sin(l, abs(m), x);
    sum = m*legendreP_over_sin + sin_times_dlegendreP;
  } else sum = sin_times_dlegendreP;

  return pref*sum;
}

Real SWSH::zeroY(const int l, const int m, const Real x)
{
  Real pref, legendreP;

  pref = pref_A_func(l, abs(m)); legendreP = plgndr(l, abs(m), x);
  if (m < 0 && m%2) pref *= -1;

  return pref*legendreP;
}

Real SWSH::dzeroYdx(const int l, const int m, const Real x)
{
  Real pref, dlegendreP;

  pref = pref_A_func(l, abs(m)); dlegendreP = dplgndr(l, abs(m), x);
  if (m < 0 && m%2) pref *= -1;

  return pref*dlegendreP;
}

Real SWSH::neg1Y(const int l, const int m, const Real x)
{
  Real pref, sin_times_dlegendreP, legendreP_over_sin;
  Real dl; // l recast as a Real - avoid integer overflows
  Real sum;

  dl = (Real)l;
  pref = pref_A_func(l, abs(m))/sqrt(dl*(dl + 1.));
  if (m < 0 && m%2) pref *= -1;
  sin_times_dlegendreP = sin_times_dplgndr(l, abs(m), x);
  if (abs(m) > 0) {
    legendreP_over_sin = plgndr_over_sin(l, abs(m), x);
    sum = m*legendreP_over_sin - sin_times_dlegendreP;
  } else sum = -sin_times_dlegendreP;

  return pref*sum;
}

Real SWSH::neg2Y(const int l, const int m, const Real x)
{
  Real pref, sinsqr_times_ddlegendreP, dlegendreP,
    legendreP_over_sinsqr;
  Real dl; // l recast as a Real - avoid integer overflows
  Real sum;

  dl = (Real)l;
  pref = pref_A_func(l, abs(m))/sqrt((dl - 1.)*dl*(dl + 1.)*(dl + 2.));
  if (m < 0 && m%2) pref *= -1;
  sinsqr_times_ddlegendreP = sinsqr_times_ddplgndr(l, abs(m), x);
  if (abs(m) > 0) {
    dlegendreP = dplgndr(l, abs(m), x);
    legendreP_over_sinsqr = plgndr_over_sinsqr(l, abs(m), x);
    sum = sinsqr_times_ddlegendreP - 2.*m*dlegendreP +
      (((Real)(m*m)) - 2.*m*x)*legendreP_over_sinsqr;
  } else sum = sinsqr_times_ddlegendreP;

  return pref*sum;
}

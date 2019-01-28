//---------------------------------------------------------------------------
//
// $Id: SWSHCGUtil.cc,v 1.9 2019/01/26 14:38:25 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// The functions in this file are used to calculate various quantities
// related to Clebsch-Gordan coefficients.  They are used by the class
// Clebsch.
//
// Scott Hughes, 24 July 1998
//

#include <cmath>
#include "Globals.h"
#include "SWSH.h"

//
// As background to these functions, let |s,l,m> represent the spherical
// (NOT spheroidal!) harmonic of spin weight s and with indices l and m.
// Further, let
//
// <s, p, m | f(\theta,\phi) | s, q, m>
//
// represent the integral
//
// $\int d\Omega {_sY^m_p}^* f(\theta,\phi) {_sY^m_q}
//

// 
// This function computes the integral
//
// <s, p, m | \cos\theta | s, q, m>
//
// See my notes on spin-weighted spheroidal harmonics for the
// formulae used to generate this.
//

Real Clebsch::xbrac(const int s, const int q, const int p, const int m)
{
  Real ans;

  ans = cgcof(q, 1, m, 0, p, m)*cgcof(q, 1, -s, 0, p, -s);
  ans *= sqrt(((Real)(2*q + 1))/((Real)(2*p + 1)));
  return(ans);
}

// 
// This function computes the integral
//
// <s,p,m | (\cos\theta)^2 | s,q,m>
//
// See my notes on spin-weighted spheroidal harmonics for the
// formulae used to generate this.
//

Real Clebsch::xsqrbrac(const int s, const int q, const int p, const int m)
{
  Real ans;

  ans = cgcof(q, 2, m, 0, p, m)*cgcof(q, 2, -s, 0, p, -s);
  ans *= (2.0/3.0)*sqrt(((Real)(2*q + 1))/((Real)(2*p + 1)));
  if (q == p) ans += (1.0/3.0);
  return(ans);
}

// 
// This function computes the integral
//
// <s,p,m | \sin\theta | s,q,mp>
//
// See Marc's and my notes for the formulae used to generate this.
//

Real Clebsch::sinthetabrac(const int s, const int p, const int q,
			   const int m, const int mp)
{
  Real ans;

  if (mp != m + 1 && mp != m - 1) {
    ans = 0.;
  } else if (mp == m + 1) {
    ans = cgcof(q, 1, mp, -1, p, m)*cgcof(q, 1, -s, 0, p, -s);
    ans *= sqrt(((Real)(4*q + 2))/((Real)(2*p + 1)));
  } else { // mp == m - 1 only possibility left ...
    ans = cgcof(q, 1, mp, 1, p, m)*cgcof(q, 1, -s, 0, p, -s);
    ans *= -sqrt(((Real)(4*q + 2))/((Real)(2*p + 1)));
  }
  return(ans);
}

//
// This function computes the Clebsch-Gordan coefficient:
// <j1, j2, m1, m2 | J, M>.  The algorithm here is from the textbook
// by Hamermesh ("Group theory and its application to physical
// problems").
//

Real Clebsch::cgcof(const int j1, const int j2, const int m1, const int m2,
		    const int J, const int M)
{
  int A, B, C, D, E;
  int lamb;
  Real rho, lnrho, CJ = 0.0, term, lnterm;

  if (m1 + m2 != M) return(0.0);
  if ((abs(j1 - j2) > J) || (j1 + j2) < J) return(0.0);
  for (lamb = 0; lamb <= (j1 + j2 - J); lamb++) {
    A = j1 + j2 - J - lamb;
    B = j1 - lamb - m1;
    C = J - j2 + lamb + m1;
    D = j2 - lamb + m2;
    E = J - j1 + lamb - m2;
    if (A < 0 || B < 0 || C < 0 || D < 0 || E < 0) 
      term = 0.0;
    else {
      lnterm = 0.5*(gsl_sf_lngamma(j1 + m1 + 1.) + gsl_sf_lngamma(j1 - m1 + 1.) +
    		    gsl_sf_lngamma(j2 + m2 + 1.) + gsl_sf_lngamma(j2 - m2 + 1.));
      lnterm -= gsl_sf_lngamma(lamb + 1.) + gsl_sf_lngamma(A + 1.) + gsl_sf_lngamma(B + 1.);
      lnterm += 0.5*(gsl_sf_lngamma(J + M + 1.) + gsl_sf_lngamma(J - M + 1.));
      lnterm -= gsl_sf_lngamma(C + 1.) + gsl_sf_lngamma(D + 1.) + gsl_sf_lngamma(E + 1.);
      term = exp(lnterm);
      if (lamb & 1) term *= -1.;
    }
    CJ += term;
  }
  lnrho = 0.5*(gsl_sf_lngamma(J - j2 + j1 + 1.) + gsl_sf_lngamma(J - j1 + j2 + 1.) +
  	       gsl_sf_lngamma(j1 + j2 - J + 1.) - gsl_sf_lngamma(j1 + j2 + J + 2.));
  rho = sqrt((Real)(2*J + 1))*exp(lnrho);
  return(rho*CJ);
}

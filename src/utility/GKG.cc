//---------------------------------------------------------------------------
//
// $Id: GKG.cc,v 1.3 2018/08/02 17:26:36 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods of class GKG.
//
// Scott Hughes, 6 January 1999.
//

#include <cmath>
#include "Globals.h"
#include "GKG.h"

#define EPSILON 1.e-9

Real GKG::RFunc(const Real r, const Real a, const Real z,
		const Real E, const Real Lz, const Real Q)
{
  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  const Real tmp2 = Kerr::Delta(r, a)*(r*r + (Lz - a*E)*(Lz - a*E) + Q);
  const Real tmp = tmp1*tmp1 - tmp2;

  return(tmp);
}

Real GKG::dr_RFunc(const Real r, const Real a, const Real z,
		   const Real E, const Real Lz, const Real Q)
{
  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  const Real tmp2 = Lz - a*E;
  const Real tmp3 = Q + r*r + tmp2*tmp2;
 
  return(4.*E*r*tmp1 - 2.*r*Kerr::Delta(r, a) - 2.*(r - 1.)*tmp3);
}

Real GKG::ddr_RFunc(const Real r, const Real a, const Real z,
		    const Real E, const Real Lz, const Real Q)
{
  const Real tmp1 = 12.*(E*E - 1.);
  const Real tmp3 = 2.*(a*a*(E*E - 1.) - Lz*Lz - Q);

  const Real ans = tmp3 + r*(12. + tmp1*r);

  return(ans);
}

Real GKG::ThetaFunc(const Real r, const Real a, const Real z,
		    const Real E, const Real Lz, const Real Q)
{
  const Real cotthsqr = z/(1. - z);
  const Real tmp = Q - cotthsqr*Lz*Lz - a*a*z*(1. - E*E);

  if(tmp < 0.) return(0.);
  else return(tmp);
}

Real GKG::PhiFunc(const Real r, const Real a, const Real z,
		  const Real E, const Real Lz, const Real Q)
{
  const Real cscthsqr = 1./(1. - z);
  const Real Delta = Kerr::Delta(r, a);
  const Real tmp = cscthsqr*Lz + a*E*((r*r + a*a)/Delta - 1.) - a*a*Lz/Delta;

  if(fabs(tmp) < EPSILON) return(0.);
  else return(tmp);
}

Real GKG::TFunc(const Real r, const Real a, const Real z,
		const Real E, const Real Lz, const Real Q)
{
  const Real Delta = Kerr::Delta(r, a);
  const Real r2pa2 = r*r + a*a;
  const Real tmp = E*((r2pa2*r2pa2)/Delta - a*a*(1. - z)) +
    a*Lz*(1. - r2pa2/Delta);

  if(fabs(tmp) < EPSILON) return(0.);
  else return(tmp);
}

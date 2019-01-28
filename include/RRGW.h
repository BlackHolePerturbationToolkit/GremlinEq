//---------------------------------------------------------------------------
//
// $Id: RRGW.h,v 1.7 2018/08/04 16:16:53 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods for calculating gravitational waves, energy and angular
// momentum losses.
//
// Scott Hughes, 17 January 1999.
//
#ifndef _RRGW_H
#define _RRGW_H

#include <cmath>
#include "Globals.h"
#include "SWSH.h"

class RRGW {
public:
  void Flux_Infinity(const Real a, const int m, const Real lamb,
		     const Real w, const Real p, const Complex ZI,
		     Real & Edot, Real & Lzdot);

  void Flux_Horizon(const Real a, const int m, const Real lamb,
		    const Real w, const Real p, const Complex ZH,
		    Real & Edot, Real & Lzdot);

  void Qdotrdot(const Real r, const Real a, const Real Q, const Real E,
		const Real Lz, const Real Edot, const Real Lzdot,
		Real & Qdot, Real & rdot);

  // Quasi-general purpose, but does not work well for inspirals!
  void Wave(const int m, const Real t_ret, const Real phi,
	    const Real S, const Real w, const Complex ZI,
	    Real & hplus, Real & hcross);

  // Good for restricted inspiral, with indices m & k.
  void Wave(const int m, const int k, const Real N_m, const Real N_k,
	    const Real phi, const Real S, const Real w,
	    const Complex ZH, Real & hplus, Real & hcross);

  // Quasi-general purpose, but does not work well for inspirals!
  void Psi4(const int m, const Real t_ret, const Real phi,
	    const Real S, const Real w, const Complex ZI,
	    Complex & psi4);

  // Used to get down-horizon fluxes.
  Real alpha_func(const Real a, const int m, const Real lamb,
		  const Real w, const Real p);
};
#endif

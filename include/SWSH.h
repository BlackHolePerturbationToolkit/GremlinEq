//---------------------------------------------------------------------------
//
// $Id: SWSH.h,v 1.15 2019/01/26 17:56:16 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes and classes I've written for the SWSH code, plus
// generically used header files.
//
// Scott Hughes, 24 July 1998
//

//
// class Clebsch contains all the routines that do things with
// Clebsch-Gordan coefficients.  See swsh_cgutil.cc for details.
//
#ifndef _SWSH_H
#define _SWSH_H

#include "Globals.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define MAXCOFS 1000
#define CHEBCOFS 256

//! Clebsch-Gordan Coefficients Class
/*! Defines methods for computing and manipulating Clebsch-Gordan coefficients */
class Clebsch {
public:
  Real xbrac(const int s, const int q, const int p, const int m);
  Real xsqrbrac(const int s, const int q, const int p, const int m);
  //
  // <s,p,m|sin(theta)|s,q,mp>
  //

  /*! Computes <s,p,m|sin(theta)|s,q,mp> */ 
  Real sinthetabrac(const int s, const int p, const int q, const int m,
		    const int mp); 

private:
  Real cgcof(const int j1, const int j2, const int m1, const int m2,
	     const int J, const int M);
};

//
// class SWSH contains all the routines that do things with
// spin-weighted spheroidal harmonics.  See swsh_spheroid.cc
// for details.
//

//! Spin-Weighted Spheroidal Harmonics Class
/*! Contains all the routines that do things with spin-weighted spheroidal harmonics */
class SWSH {
public:
  //! A constructor for the SWSH Class.
  /*!
      Initialized with (s,l,m,spheroidicity)
  */
  SWSH(const int ss, const int ll, const int mm, const Real spintimesfreq);

  void expand(Real *E, Real *b, int *n);
  Real error(const Real E, const Real b[], const int n);

  Real spheroid(const Real x);
  Real l2dagspheroid(const Real x);
  Real l1dagl2dagspheroid(const Real x);
  Complex edthbaredthbarspheroid(const Real a, const Real x);

  Real pos2Y(const int l, const int m, const Real x);
  Real pos1Y(const int l, const int m, const Real x);
  Real zeroY(const int l, const int m, const Real x);
  Real dzeroYdx(const int l, const int m, const Real x);
  Real neg1Y(const int l, const int m, const Real x);
  Real neg2Y(const int l, const int m, const Real x);

  Real E, lambda;

  Real b[MAXCOFS + 1];
  int lmin; 						//!< Minimum {\it harmonic} index in array
  int N; 							//!< Maximum index in array

  int Gl;

private:
  Real pref_A_func(const int l, const int m);

  Real plgndr(const int l, const int m, const Real x);
  Real plgndr_over_sin(const int l, const int m, const Real x);
  Real plgndr_over_sinsqr(const int l, const int m, const Real x);

  Real dplgndr(const int l, const int m, const Real x);
  Real sin_times_dplgndr(const int l, const int m, const Real x);

  Real sinsqr_times_ddplgndr(const int l, const int m, const Real x);
  //
  // "Raw" spheroidal routines.
  //
  Real raw_spheroid(const Real x);
  Real raw_l2dagspheroid(const Real x);
  Real raw_l1dagl2dagspheroid(const Real x);
  //
  // Routines to compute Chebyshev expansion coefficients.
  //
  void cheb_spheroid();
  void cheb_l2dagspheroid();
  void cheb_l1dagl2dagspheroid();
  //
  // Storage of Chebyshev expansion coefficients.
  //
  Real spheroid_cofs[CHEBCOFS];
  Real l2dagspheroid_cofs[CHEBCOFS];
  Real l1dagl2dagspheroid_cofs[CHEBCOFS];
  //
  // Last nontrivial expansion coefficients.
  //
  int N_spheroid;
  int N_l2dagspheroid;
  int N_l1dagl2dagspheroid;

  int findeig(Real *M, const int K, Real *e, Real b[]);
  void load_M(Real *M, const int K);

  // Used globally by the functions:
  //
  // s = spin weight; Gl,Gm = the spheroidal indices (attributed
  // globally); aw = a*omega; lmin = min(s,Gm), the minimum allowed l;
  // *b will be allocated to an array of expansion coefficients;
  // E = spheroidal eigenvalue; N = number of coefficients in
  // expansion.
  int s, Gm;
  Real aw;
};

#endif

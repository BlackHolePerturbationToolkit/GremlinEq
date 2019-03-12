//---------------------------------------------------------------------------
//
// $Id: TidalH.h,v 1.13 2019/01/27 09:43:08 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// All the stuff we need for doing analyses of tidally distorted
// horizons, including embeddings.
//
// Scott A. Hughes, 30 May 2014 (started)
//
#ifndef _TIDALHORIZ_H
#define _TIDALHORIZ_H

#include "Globals.h"
#include "SWSH.h"

#define mRq_IGRND 1
#define mD_C_IGRND 2
#define mD_D_IGRND 3
//! Tidal Horizon Class
/*! This class defines functions which are useful for analyzing how the on-horizon Teukolsky solution affects the geometry of a black hole. This methods were extensively used in work with former MIT student Stephen O’Sullivan, but may not be of broad interest; as such, they may be removed from general release. */
class TidalH {
public:
  TidalH(const Real spin, const int ellmax);
  ~TidalH();
  Real a, rp, eps, Kph;
  int lmax;
  SWSH **swshp;
  Complex **ZH, **Clm, **Epslm;
  Real *w, *p, **lamb;

  // Global variables k, l, m
  int Gq, Gl, Gm;

  int INTEGRAND;

  void loadClm();
  void loadEpslm();
  Complex mR_vector(const int q, const int m);
  Real mD_matrix(const int q, const int ell, const int m);

  Real R1lm(const int l, const int m, const Real x, const Real psi);
  Real epsr(const int l, const int m, const Real x, const Real psi);

  Real C0(const Real x), C1(const Real x);
  Real D(const Real x);

  Complex IGRND(const Real x);
  //
/*   Real IGRNDre(Real x); */
  static Real IGRNDre_wrapper(Real x, void *params) {
    return static_cast<TidalH*>(params)->IGRND(x).real();
  };
  //
/*   Real IGRNDim(Real x); */
  static Real IGRNDim_wrapper(Real x, void *params) {
    return static_cast<TidalH*>(params)->IGRND(x).imag();
  };
  //
};
#endif

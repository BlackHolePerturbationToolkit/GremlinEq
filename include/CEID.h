//---------------------------------------------------------------------------
//
// $Id: CEID.h,v 1.33 2019/01/24 20:47:07 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial inspiral data.
//
#ifndef _CEID_H
#define _CEID_H

#include <cmath>
#include "Globals.h"
#include "CEDR.h"
#include "SWSH.h"
#include "RRGW.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define LRANGE 100
//! Circular Equatorial Inspiral Data Class
/*! This class uses CEDR, and provides methods for smoothly interpolating through a large set of output from many orbits to study waveforms and inspirals. */
class CEID {
public:
  CEID(char *basename, const Real arrmin, const Real rmax, const Real dr,
       const int ellmax);
  CEID(char *basename, const int Nhi, const int ellmax);
  ~CEID();

  Real rdot, rdot_noH, Omega;
  Real Omega_max, rmin;
  Real EdotI, LzdotI;
  Real EdotH, LzdotH;
  Real hp, hc;
  Real r_isco, a;
  //
  // Loads spheroid data.
  //
  void Spheroids(const Real costheta_view);
  //
  // Get stuff, interpolated to current coordinate in parameter space
  //
  void Get_rdot_omega(const Real r);
  void Get_fluxes(const Real r);
  void Get_wave(const Real r, const Real phi, const Real phi_view);
  void Get_flux_Smode(const Real r, const int l, const int m,
		      Real & EdotHlm, Real & EdotIlm);
  void Get_flux_mmode(const Real r, const int m,
		      Real & EdotHm, Real & EdotIm);
  void Get_flux_lmode(const Real r, const int l,
		      Real & EdotHl, Real & EdotIl);
  void Get_flux_Ymode(const Real r, const int l, const int m,
		      Real & EdotHlm, Real & EdotIlm);
  void Get_ZIlm(const Real r, const int l, const int m, Complex & ZIlm);
  void Get_CIlm(const Real r, const Real phi, const int l, const int m,
		Complex & CIlm);
  //
  // Before using Get_CIlm or Get_flux_Ymode, you must run Spline_CIlm()!
  // Would be best to put this in the constructor, but it is so CPU
  // expensive we leave it out to be run for a particular mode.
  //
  void Spline_CIlm(const int l, const int m);

  int jmax;
  Real *r_arr;

  int *max_l_computed; // The largest value of l computed for a given v.

private:
  void Allocate();
  void ReadIn(char *inname);
  void Make_radial_splines();

  int j, lmax, l, m;
  //
  Real *rstar_arr;
  //
  Real *Omega_arr, *rdot_arr, *rdot_noH_arr;
  gsl_interp_accel *Omega_acc, *rdot_acc, *rdot_noH_acc;
  gsl_spline *Omega_spl, *rdot_spl, *rdot_noH_spl;
  //
  Real *EdotI_arr, *LzdotI_arr, *EdotH_arr, *LzdotH_arr;
  gsl_interp_accel *EdotI_acc, *LzdotI_acc, *EdotH_acc, *LzdotH_acc;
  gsl_spline *EdotI_spl, *LzdotI_spl, *EdotH_spl, *LzdotH_spl;
  //
  // Used for waveforms.
  //
  Real ***ZI_re, ***ZI_im, ***CI_re, ***CI_im, ***ZH_re, ***ZH_im, ***Spheroid;
  gsl_interp_accel ***ZI_re_acc, ***ZI_im_acc, ***CI_re_acc, ***CI_im_acc;
  gsl_interp_accel ***ZH_re_acc, ***ZH_im_acc, ***Spheroid_acc;
  gsl_spline ***ZI_re_spl, ***ZI_im_spl, ***CI_re_spl, ***CI_im_spl;
  gsl_spline ***ZH_re_spl, ***ZH_im_spl, ***Spheroid_spl;
};
#endif

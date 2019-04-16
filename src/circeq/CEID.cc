//---------------------------------------------------------------------------
//
// $Id: CEID.cc,v 1.46 2019/01/24 20:41:23 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial inspiral data.
//
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "Tensors.h"
#include "CEDR.h"
#include "CEID.h"
#include "SWSH.h"
#include "RRGW.h"

CEID::CEID(char *basename, const Real arrmin, const Real rmax, const Real dr,
	   const int ellmax) : rmin(arrmin), lmax(ellmax)
{
  char inname[150];
  //
  // This assignment of jmax is a bit whacky, but it works well.
  //
  jmax = (int)ceil((rmax - rmin)/dr + 0.0001);
  //
  // Allocate memory
  //
  Allocate();
  //
  // Read in data.
  //
  Omega_max = 0.;
  for (j = 1; j <= jmax; j++) {
    Real radius = rmin + (j - 1)*dr;
    if (dr < 0.01)
      sprintf(inname, "%s_r%4.3lf.h5", basename, radius);
    else if (dr < 0.1)
      sprintf(inname, "%s_r%3.2lf.h5", basename, radius);
    else
      sprintf(inname, "%s_r%2.1lf.h5", basename, radius);
    ReadIn(inname);
  }
  //
  // Make radial splines
  //
  Make_radial_splines();
}

CEID::CEID(char *basename, const int Nhi, const int ellmax) : jmax(Nhi), lmax(ellmax)
{
  char inname[150];
  //
  // Allocate memory
  //
  Allocate();
  //
  // Read in data.
  //
  Omega_max = 0.;
  //
  // FIX THIS: Internal data arrays are zero offset; external data is unit offset.
  //
  for (j = 0; j < jmax; j++) {
    sprintf(inname, "%s_%d.h5", basename, j+1);
    ReadIn(inname);
  }
  //
  // Make radial splines
  //
  Make_radial_splines();
  rmin = r_arr[1];
}

void CEID::Allocate()
{
  //
  // The stuff we need for the gsl spline routines
  //
  r_arr = Tensor<Real>::vector(0, jmax-1);
  rstar_arr = Tensor<Real>::vector(0, jmax-1);
  Omega_arr = Tensor<Real>::vector(0, jmax-1);
  rdot_arr = Tensor<Real>::vector(0, jmax-1);
  rdot_noH_arr = Tensor<Real>::vector(0, jmax-1);
  EdotI_arr = Tensor<Real>::vector(0, jmax-1);
  LzdotI_arr = Tensor<Real>::vector(0, jmax-1);
  EdotH_arr = Tensor<Real>::vector(0, jmax-1);
  LzdotH_arr = Tensor<Real>::vector(0, jmax-1);
  //
  Omega_acc = gsl_interp_accel_alloc();
  rdot_acc = gsl_interp_accel_alloc();
  rdot_noH_acc = gsl_interp_accel_alloc();
  Omega_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  rdot_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  rdot_noH_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  //
  EdotI_acc = gsl_interp_accel_alloc();
  EdotH_acc = gsl_interp_accel_alloc();
  LzdotI_acc = gsl_interp_accel_alloc();
  LzdotH_acc = gsl_interp_accel_alloc();
  EdotI_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  EdotH_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  LzdotI_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  LzdotH_spl = gsl_spline_alloc(gsl_interp_cspline, jmax);
  //
  max_l_computed = Tensor<int>::vector(0, jmax-1);
  //
  if (lmax > 0) {
    ZI_re = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    ZI_im = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    CI_re = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    CI_im = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    ZH_re = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    ZH_im = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    Spheroid = Tensor<Real>::tensor3(2, lmax, -lmax, lmax, 0, jmax-1);
    //
    ZI_re_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    ZI_re_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    ZI_im_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    ZI_im_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    CI_re_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    CI_re_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    CI_im_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    CI_im_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    ZH_re_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    ZH_re_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    ZH_im_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    ZH_im_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    Spheroid_acc = Tensor<gsl_interp_accel>::matrixptr(2, lmax, -lmax, lmax);
    Spheroid_spl = Tensor<gsl_spline>::matrixptr(2, lmax, -lmax, lmax);
    //
    // Finish allocating memory to the splines; preload the data with
    // zeros
    //
    for (l = 2; l <= lmax; l++) {
      for (m = -l; m <= l; m++) {
 	ZI_re_acc[l][m] = gsl_interp_accel_alloc();
 	ZI_re_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	ZI_im_acc[l][m] = gsl_interp_accel_alloc();
 	ZI_im_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	CI_re_acc[l][m] = gsl_interp_accel_alloc();
 	CI_re_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	CI_im_acc[l][m] = gsl_interp_accel_alloc();
 	CI_im_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	ZH_re_acc[l][m] = gsl_interp_accel_alloc();
 	ZH_re_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	ZH_im_acc[l][m] = gsl_interp_accel_alloc();
 	ZH_im_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
 	Spheroid_acc[l][m] = gsl_interp_accel_alloc();
 	Spheroid_spl[l][m] = gsl_spline_alloc(gsl_interp_cspline, jmax);
	for(j = 0; j < jmax; j++) {
	  ZI_re[l][m][j] = 0.;
	  ZI_im[l][m][j] = 0.;
	  CI_re[l][m][j] = 0.;
	  CI_im[l][m][j] = 0.;
	  ZH_re[l][m][j] = 0.;
	  ZH_im[l][m][j] = 0.;
	  Spheroid[l][m][j] = 0.;
	}
      }
    }
  }
}

void CEID::ReadIn(char *inname)
{
  CEDR cedr(inname);
  r_arr[j] = cedr.r;
  a = cedr.a;
  Real Lz = cedr.Lz;
  rstar_arr[j] = Kerr::rstar(r_arr[j], a);
  Omega_arr[j] = cedr.Om_phi;
  if (fabs(Omega_arr[j]) > Omega_max)
    Omega_max = fabs(Omega_arr[j]);
  //
  Real EdotI_in = 0., EdotH_in = 0.;
  Real LzdotI_in = 0., LzdotH_in = 0.;
  Real rdotI_in = 0., rdotH_in = 0.;
  //
  // Read in data for each mode
  //
  max_l_computed[j] = 0;
  for (l = cedr.lmin; l <= cedr.lmax; l++) {
    if (l > max_l_computed[j]) max_l_computed[j] = l;
    for (m = cedr.mmin; m <= cedr.mmax; m++) {
      if (cedr.ReadData(l, m)) {
	//
	// Factor of 2 because data only includes positive m, but we've
	// got symmetry.  Note that Edot and Lzdot are the fluxes
	// carried by the radiation, NOT the change in the particle's
	// energy and angular momentum.
	//
	EdotI_in += 2.*cedr.EdotI;
	EdotH_in += 2.*cedr.EdotH;
	LzdotI_in += 2.*cedr.LzdotI;
	LzdotH_in += 2.*cedr.LzdotH;
	rdotI_in += 2.*cedr.rdotI;
	rdotH_in += 2.*cedr.rdotH;
	if (l <= lmax) {
	  ZI_re[l][m][j] = cedr.ZI.real();
	  ZI_im[l][m][j] = cedr.ZI.imag();
	  ZH_re[l][m][j] = cedr.ZH.real();
	  ZH_im[l][m][j] = cedr.ZH.imag();
	  Real sign = (l%2) ? -1. : 1.;
	  ZI_re[l][-m][j] = sign*ZI_re[l][m][j];
	  ZI_im[l][-m][j] = -sign*ZI_im[l][m][j];
	  ZH_re[l][-m][j] = sign*ZH_re[l][m][j];
	  ZH_im[l][-m][j] = -sign*ZH_im[l][m][j];
	} // if (l
      } // if (cedr.ReadData
    } // for (m
  } // for (l
  EdotI_arr[j] = EdotI_in; EdotH_arr[j] = EdotH_in;
  LzdotI_arr[j] = LzdotI_in; LzdotH_arr[j] = LzdotH_in;
  //
  // Use r_isco to clear out divergence in rdot towards the lso
  //
  if (Lz > 0.) { // prograde
    r_isco = Kerr::isco_pro(a);
  } else { // retrograde
    r_isco = Kerr::isco_ret(a);
  }
  rdot_arr[j] = (rdotI_in + rdotH_in)*(r_arr[j] - r_isco);
  rdot_noH_arr[j] = rdotI_in*(r_arr[j] - r_isco);
}

void CEID::Make_radial_splines()
{
  gsl_spline_init(Omega_spl, rstar_arr, Omega_arr, jmax);
  gsl_spline_init(rdot_spl, rstar_arr, rdot_arr, jmax);
  gsl_spline_init(rdot_noH_spl, rstar_arr, rdot_noH_arr, jmax);
  gsl_spline_init(EdotI_spl, rstar_arr, EdotI_arr, jmax);
  gsl_spline_init(EdotH_spl, rstar_arr, EdotH_arr, jmax);
  gsl_spline_init(LzdotI_spl, rstar_arr, LzdotI_arr, jmax);
  gsl_spline_init(LzdotH_spl, rstar_arr, LzdotH_arr, jmax);
  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      gsl_spline_init(ZI_re_spl[l][m], rstar_arr, ZI_re[l][m], jmax);
      gsl_spline_init(ZI_im_spl[l][m], rstar_arr, ZI_im[l][m], jmax);
      gsl_spline_init(ZH_re_spl[l][m], rstar_arr, ZH_re[l][m], jmax);
      gsl_spline_init(ZH_im_spl[l][m], rstar_arr, ZH_im[l][m], jmax);
    }
  }
}

void CEID::Spline_CIlm(const int l, const int m)
{
  int k, j, lmin, ell;
  //
  // First, make an array of CIlm for each location on our grid.
  //
  for (j = 0; j < jmax; j++) {
    Complex CIlm_here = Complex(0., 0.);
    Real where = Omega_arr[j];
    lmin = Max(2, abs(m));
    //
    k = l;
    for (ell = lmin; ell <= lmax; ell++) {
      SWSH swsh_neg2_ell(-2, ell, m, a*((Real)m)*where);
      if (k - lmin + 1 <= swsh_neg2_ell.N)
	CIlm_here += (ZI_re[ell][m][j] + II*ZI_im[ell][m][j])
	  *swsh_neg2_ell.b[k - lmin + 1];
    }
    CI_re[l][m][j] = CIlm_here.real();
    CI_im[l][m][j] = CIlm_here.imag();
  }
  gsl_spline_init(CI_re_spl[l][m], rstar_arr, CI_re[l][m], jmax);
  gsl_spline_init(CI_im_spl[l][m], rstar_arr, CI_im[l][m], jmax);
}

CEID::~CEID()
{
  Tensor<int>::free_vector(max_l_computed, 0, jmax-1);
  //
  Tensor<Real>::free_vector(r_arr, 0, jmax-1);
  Tensor<Real>::free_vector(rstar_arr, 0, jmax-1);
  Tensor<Real>::free_vector(Omega_arr, 0, jmax-1);
  Tensor<Real>::free_vector(rdot_arr, 0, jmax-1);
  Tensor<Real>::free_vector(rdot_noH_arr, 0, jmax-1);
  Tensor<Real>::free_vector(EdotI_arr, 0, jmax-1);
  Tensor<Real>::free_vector(EdotH_arr, 0, jmax-1);
  Tensor<Real>::free_vector(LzdotI_arr, 0, jmax-1);
  Tensor<Real>::free_vector(LzdotH_arr, 0, jmax-1);
  //
  gsl_spline_free(Omega_spl);
  gsl_spline_free(rdot_spl);
  gsl_spline_free(rdot_noH_spl);
  gsl_interp_accel_free(Omega_acc);
  gsl_interp_accel_free(rdot_acc);
  gsl_interp_accel_free(rdot_noH_acc);
  gsl_spline_free(EdotI_spl);
  gsl_spline_free(EdotH_spl);
  gsl_spline_free(LzdotI_spl);
  gsl_spline_free(LzdotH_spl);
  gsl_interp_accel_free(EdotI_acc);
  gsl_interp_accel_free(EdotH_acc);
  gsl_interp_accel_free(LzdotI_acc);
  gsl_interp_accel_free(LzdotH_acc);
  if (lmax > 0) {
    for (l = 2; l <= lmax; l++) {
      for (m = -lmax; m <= lmax; m++) {
	gsl_interp_accel_free(ZI_re_acc[l][m]);
	gsl_interp_accel_free(ZI_im_acc[l][m]);
	gsl_spline_free(ZI_re_spl[l][m]);
	gsl_spline_free(ZI_im_spl[l][m]);
	gsl_interp_accel_free(CI_re_acc[l][m]);
	gsl_interp_accel_free(CI_im_acc[l][m]);
	gsl_spline_free(CI_re_spl[l][m]);
	gsl_spline_free(CI_im_spl[l][m]);
	gsl_interp_accel_free(ZH_re_acc[l][m]);
	gsl_interp_accel_free(ZH_im_acc[l][m]);
	gsl_spline_free(ZH_re_spl[l][m]);
	gsl_spline_free(ZH_im_spl[l][m]);
	gsl_interp_accel_free(Spheroid_acc[l][m]);
	gsl_spline_free(Spheroid_spl[l][m]);
      }
    }
    Tensor<Real>::free_tensor3(ZI_re, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<Real>::free_tensor3(ZI_im, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<gsl_interp_accel>::free_matrixptr(ZI_re_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_interp_accel>::free_matrixptr(ZI_im_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(ZI_re_spl, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(ZI_im_spl, 2, lmax, -lmax, lmax);
    Tensor<Real>::free_tensor3(CI_re, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<Real>::free_tensor3(CI_im, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<gsl_interp_accel>::free_matrixptr(CI_re_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_interp_accel>::free_matrixptr(CI_im_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(CI_re_spl, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(CI_im_spl, 2, lmax, -lmax, lmax);
    Tensor<Real>::free_tensor3(ZH_re, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<Real>::free_tensor3(ZH_im, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<gsl_interp_accel>::free_matrixptr(ZH_re_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_interp_accel>::free_matrixptr(ZH_im_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(ZH_re_spl, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(ZH_im_spl, 2, lmax, -lmax, lmax);
    Tensor<Real>::free_tensor3(Spheroid, 2, lmax, -lmax, lmax, 0, jmax-1);
    Tensor<gsl_interp_accel>::free_matrixptr(Spheroid_acc, 2, lmax, -lmax, lmax);
    Tensor<gsl_spline>::free_matrixptr(Spheroid_spl, 2, lmax, -lmax, lmax);
  }
}

void CEID::Spheroids(const Real costheta_view)
{
  if (lmax > 0) {
    for (j = 0; j < jmax; j++) {
      for (l = 2; l <= lmax; l++) {
	for (m = -l; m <= l; m++) {
	  SWSH swsh(-2, l, m, a*m*Omega_arr[j]);
	  Spheroid[l][m][j] = swsh.spheroid(costheta_view);
	}
      }
    }
    for (l = 2; l <= lmax; l++)
      for (m = -l; m <= l; m++)
	gsl_spline_init(Spheroid_spl[l][m], rstar_arr, Spheroid[l][m], jmax);
  } else return;
}

void CEID::Get_rdot_omega(const Real r)
{
  const Real rs = Kerr::rstar(r, a);

  rdot = gsl_spline_eval(rdot_spl, rs, rdot_acc);
  rdot /= (r - r_isco);
  rdot_noH = gsl_spline_eval(rdot_noH_spl, rs, rdot_noH_acc);
  rdot_noH /= (r - r_isco);
  Omega = gsl_spline_eval(Omega_spl, rs, Omega_acc);
}

void CEID::Get_fluxes(const Real r)
{
  const Real rs = Kerr::rstar(r, a);

  EdotI = gsl_spline_eval(EdotI_spl, rs, EdotI_acc);
  EdotH = gsl_spline_eval(EdotH_spl, rs, EdotH_acc);
  LzdotI = gsl_spline_eval(LzdotI_spl, rs, LzdotI_acc);
  LzdotH = gsl_spline_eval(LzdotH_spl, rs, LzdotH_acc);
}

void CEID::Get_wave(const Real r, const Real phi, const Real phi_view)
{
  const Real rs = Kerr::rstar(r, a);

  hp = 0.; hc = 0.;
  Real ZI_re_int, ZI_im_int, Spheroid_int, dhc, dhp;
  Complex Zhere;
  RRGW rrgw;
  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      if (m != 0) {
	ZI_re_int = gsl_spline_eval(ZI_re_spl[l][m], rs, ZI_re_acc[l][m]);
	ZI_im_int = gsl_spline_eval(ZI_im_spl[l][m], rs, ZI_im_acc[l][m]);
	Spheroid_int = gsl_spline_eval(Spheroid_spl[l][m], rs, Spheroid_acc[l][m]);
	Zhere = ZI_re_int + II*ZI_im_int;
	//
	// This goes to the second verson of "Wave" in rrgw --- arguments
	// are m, k, number of axial cycles, number of polar cycles, viewing
	// angle phi, spheroid, omega_{mk}, ZI_{lmk}, dhplus, dhcross.
	// The last two arguments are returned non-zero.
	rrgw.Wave(m, 0, phi/(2.*M_PI), 0., phi_view, Spheroid_int, m*Omega, Zhere, dhp, dhc);
	hp += dhp; hc += dhc;
      }
    }
  }
}

//
// Provides the flux in the (l, m) spheroidal mode.
void CEID::Get_flux_Smode(const Real r, const int l, const int m,
			  Real & EdotHlm, Real & EdotIlm)
{
  const Real rs = Kerr::rstar(r, a);
  Complex ZI, ZH;
  Real ZI_re_int, ZI_im_int;
  Real ZH_re_int, ZH_im_int;
  Real junk;
  RRGW rrgw;

  ZI_re_int = gsl_spline_eval(ZI_re_spl[l][m], rs, ZI_re_acc[l][m]);
  ZI_im_int = gsl_spline_eval(ZI_im_spl[l][m], rs, ZI_im_acc[l][m]);
  ZI = ZI_re_int + II*ZI_im_int;
  ZH_re_int = gsl_spline_eval(ZH_re_spl[l][m], rs, ZH_re_acc[l][m]);
  ZH_im_int = gsl_spline_eval(ZH_im_spl[l][m], rs, ZH_im_acc[l][m]);
  ZH = ZH_re_int + II*ZH_im_int;
  //
  // A few useful global defines.
  Real w = ((Real)m)*Omega;
  Real p = w - 0.5*m*a/Kerr::rplus(a);
  SWSH swsh(-2, l, m, a*w);
  //
  rrgw.Flux_Infinity(a, m, swsh.lambda, w, p, ZI, EdotIlm, junk);
  rrgw.Flux_Horizon(a, m, swsh.lambda, w, p, ZH, EdotHlm, junk);
}

//
// Provides the flux in the m azimuthal mode.
void CEID::Get_flux_mmode(const Real r, const int m,
			  Real & EdotHm, Real & EdotIm)
{
  const Real rs = Kerr::rstar(r, a);
  int l;
  Complex ZI, ZH;
  Real ZI_re_int, ZI_im_int;
  Real ZH_re_int, ZH_im_int;
  Real EdotHlm, EdotIlm;
  Real junk;
  RRGW rrgw;

  EdotHm = 0.;
  EdotIm = 0.;
  int lmin = Max(2, abs(m));
  for (l = lmin; l <= lmax; l++) {
    ZI_re_int = gsl_spline_eval(ZI_re_spl[l][m], rs, ZI_re_acc[l][m]);
    ZI_im_int = gsl_spline_eval(ZI_im_spl[l][m], rs, ZI_im_acc[l][m]);
    ZI = ZI_re_int + II*ZI_im_int;
    ZH_re_int = gsl_spline_eval(ZH_re_spl[l][m], rs, ZH_re_acc[l][m]);
    ZH_im_int = gsl_spline_eval(ZH_im_spl[l][m], rs, ZH_im_acc[l][m]);
    ZH = ZH_re_int + II*ZH_im_int;
    //
    Real w = ((Real)m)*Omega;
    Real p = w - 0.5*m*a/Kerr::rplus(a);
    SWSH swsh(-2, l, m, a*w);
    //
    rrgw.Flux_Infinity(a, m, swsh.lambda, w, p, ZI, EdotIlm, junk);
    rrgw.Flux_Horizon(a, m, swsh.lambda, w, p, ZH, EdotHlm, junk);
    EdotHm += EdotHlm;
    EdotIm += EdotIlm;
  }
}

//
// Provides the flux in the l spheroidal mode
void CEID::Get_flux_lmode(const Real r, const int l,
			  Real & EdotHl, Real & EdotIl)
{
  const Real rs = Kerr::rstar(r, a);
  int m;
  Complex ZI, ZH;
  Real ZI_re_int, ZI_im_int;
  Real ZH_re_int, ZH_im_int;
  Real EdotHlm, EdotIlm;
  Real junk;
  RRGW rrgw;

  EdotHl = 0.;
  EdotIl = 0.;
  for (m = 1; m <= l; m++) {
    ZI_re_int = gsl_spline_eval(ZI_re_spl[l][m], rs, ZI_re_acc[l][m]);
    ZI_im_int = gsl_spline_eval(ZI_im_spl[l][m], rs, ZI_im_acc[l][m]);
    ZI = ZI_re_int + II*ZI_im_int;
    ZH_re_int = gsl_spline_eval(ZH_re_spl[l][m], rs, ZH_re_acc[l][m]);
    ZH_im_int = gsl_spline_eval(ZH_im_spl[l][m], rs, ZH_im_acc[l][m]);
    ZH = ZH_re_int + II*ZH_im_int;
    //
    Real w = ((Real)m)*Omega;
    Real p = w - 0.5*m*a/Kerr::rplus(a);
    SWSH swsh(-2, l, m, a*w);
    //
    rrgw.Flux_Infinity(a, m, swsh.lambda, w, p, ZI, EdotIlm, junk);
    rrgw.Flux_Horizon(a, m, swsh.lambda, w, p, ZH, EdotHlm, junk);
    EdotHl += 2.*EdotHlm;
    EdotIl += 2.*EdotIlm;
  }
}

//
// Provides the flux in the (l, m) spherical mode.
void CEID::Get_flux_Ymode(const Real r, const int l, const int m,
			  Real & EdotHlm, Real & EdotIlm)
{
  const Real rs = Kerr::rstar(r, a);
  Complex CI, *ZH;
  Real CI_re_int, CI_im_int;
  Real ZH_re_int, ZH_im_int;
  int j, k, ell, lmin;
  Real junk;
  RRGW rrgw;

  CI_re_int = gsl_spline_eval(CI_re_spl[l][m], rs, CI_re_acc[l][m]);
  CI_im_int = gsl_spline_eval(CI_im_spl[l][m], rs, CI_im_acc[l][m]);
  CI = CI_re_int + II*CI_im_int;
  //
  // A few useful global defines.
  Real w = ((Real)m)*Omega;
  Real p = w - 0.5*m*a/Kerr::rplus(a);
  Real rp = Kerr::rplus(a);
  Real eps = sqrt((1. - a)*(1. + a))/(4.*rp);
  //
  EdotIlm = abs(CI)*abs(CI)/(4.*M_PI*w*w);
  //
  // Now, the down-horizon stuff.  This is a bit complicated, and
  // does not correspond nicely to the methods in RRGW.
  //
  ZH = Tensor<Complex>::vector(2, lmax);
  // Index m is implied, since everything is done here at the same m.
  for (k = 2; k <= lmax; k++) {
    ZH_re_int = gsl_spline_eval(ZH_re_spl[l][m], rs, ZH_re_acc[l][m]);
    ZH_im_int = gsl_spline_eval(ZH_im_spl[l][m], rs, ZH_im_acc[l][m]);
    ZH[k] = ZH_re_int + II*ZH_im_int;
  }
  Complex UHlm = Complex(0., 0.);
  lmin = Max(2, abs(m));
  //
  j = l; // Just for correspondance with my notes.
  for (ell = lmin; ell <= lmax; ell++) {
    SWSH swsh_neg2_ell(-2, ell, m, a*w);
    SWSH swsh_pos2_ell(2, ell, m, a*w);
    if (j - lmin + 1 <= swsh_pos2_ell.N) {
      const Real lamb = swsh_neg2_ell.lambda;
      const Real lp2 = lamb + 2.;
      
      const Real aw_m_m = a*w - ((Real)m);
      const Real tmp1 = lp2*lp2 - 4.*a*w*aw_m_m;
      const Real tmp2 = lamb*lamb - 36.*a*w*aw_m_m;
      const Real tmp3 = 48.*a*w*(2.*lamb + 3.)*(2.*a*w - ((Real)m));
      
      const Real abs_sqr_clm = tmp1*tmp2 + tmp3 + 144.*w*w*(1. - a*a);
      const Real Im_clm = 12.*w;
      const Real Re_clm = sqrt(abs_sqr_clm - Im_clm*Im_clm);
      const Complex clm = Complex(Re_clm, Im_clm);
      
      const Complex beta_numer = 64.*16.*pow(Kerr::rplus(a), 4.)*II*p*
	(p*p + 4.*eps*eps)*(4.*eps - II*p);
      const Complex beta = beta_numer/clm;
      UHlm += beta*ZH[ell]*swsh_pos2_ell.b[j - lmin + 1];
    }
  }
  EdotHlm = w*abs(UHlm)*abs(UHlm)/
    (64.*M_PI*p*(p*p + 4.*eps*eps)*pow(2.*rp, 3.));
  Tensor<Complex>::free_vector(ZH, 2, lmax);
}

void CEID::Get_ZIlm(const Real r, const int l, const int m, Complex & ZIlm)
{
  //
  // Assumes that Spline_Clm() has been run for these values of l!
  //
  const Real rs = Kerr::rstar(r, a);
  Real ZI_re_int, ZI_im_int;

  ZI_re_int = gsl_spline_eval(ZI_re_spl[l][m], rs, ZI_re_acc[l][m]);
  ZI_im_int = gsl_spline_eval(ZI_im_spl[l][m], rs, ZI_im_acc[l][m]);

  ZIlm = ZI_re_int + II*ZI_im_int;
}

void CEID::Get_CIlm(const Real r, const Real phi, const int l, const int m,
		    Complex & CIlm)
{
  //
  // Assumes that Spline_Clm() has been run for these values of l!
  //
  const Real rs = Kerr::rstar(r, a);
  Real CI_re_int, CI_im_int;

  CI_re_int = gsl_spline_eval(CI_re_spl[l][m], rs, CI_re_acc[l][m]);
  CI_im_int = gsl_spline_eval(CI_im_spl[l][m], rs, CI_im_acc[l][m]);

  CIlm = (CI_re_int + II*CI_im_int)*exp(-II*((Real)m)*phi);
}

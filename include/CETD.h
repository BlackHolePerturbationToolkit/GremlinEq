//---------------------------------------------------------------------------
//
// $Id: CETD.h,v 1.13 2018/08/04 15:56:42 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial Teukolsky driver
// Scott Hughes, 26 July 2003
//
#ifndef _CETD_H
#define _CETD_H

#include <cmath>
#include <hdf5_hl.h>
#include "Globals.h"
#include "CEKG.h"
#include "SWSH.h"
#include "FT.h"
#include "CEKR.h"
#include "RRGW.h"

//! Circular Equatorial Teukolsky Driver Class
/*! A class which is used to “drive” studies that look at a many circular equatorial orbits. 
It is essentially a holder for methods that solve CEKR repeatedly. This class is also used to write the HDF5 data files (the “d” stands for data as well as driver). 
The synopsis is that the HDF5 file contains 2 groups: the 2 modes (i.e., the indices l and m) and the parameters (orbital radius r, spin parameter a, energy E, axial angular momentum Lz, and axial frequency Om_phi). 
Into the modes group goes 19 different bits of data,
* - The spin-weighted spheroidal harmonic eigenvalue lambda
* - The real and the imaginary parts of Z_Inf (the amplitude of the “to infinity” Teukolsky amplitude)
* - The real and the imaginary parts of Z_H (the amplitude of the “down horizon” Teukolsky amplitude)
* - The real and the imaginary parts of Rin (the ingoing separated radial part of the Teukolsky function, at the orbit)
* - The real and the imaginary parts of d/dr(Rin)
* - The real and the imaginary parts of Rup (the outgoing separated radial part of the Teukolsky function, at the orbit)
* - The real and the imaginary parts of d/dr(Rup)
* - 4 fluxes carried by the radiation, Edot_Inf, Edot_H, Lzdot_Inf, Lzdot_H
* - The rate of change of orbital radius arising from radiation flux to infinity, rdot_Inf
* - The rate of change of orbital radius arising from the down-horizon radiation flux, rdot_H
*/
class CETD {
public:
  CETD(const int orbitsense, const Real rad, const Real spin, char outbase[]);
  ~CETD();
  //
  void Driver(const int lmax);
  void Driver(const Real EPS_L, const int lmax_min);
  void Driver(const Real EPS_L, const int lmax, const int lmax_min);
  //
  void DoHarmonic(const int l, const int m);
  //
  char outname[256];
  hid_t hdffile;
  //
  RRGW *rrgw;
  CEKG *cekg;
  //
  int proret;
  int l, m;
  Real r, a;
  //
  // Range of harmonics computed by the driver.
  //
  int lmin, lmax, mmin, mmax;
  //
  Real EdotH, EdotI;
  Complex ZH, ZI;
};
#endif

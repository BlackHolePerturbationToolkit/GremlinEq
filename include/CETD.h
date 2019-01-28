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

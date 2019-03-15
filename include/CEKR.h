//---------------------------------------------------------------------------
//
// $Id: CEKR.h,v 1.4 2018/08/02 23:04:34 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Class CEKR (Circular, Equatorial Kerr Radiation).
//
// Scott Hughes, 11 January 1999.
//
#ifndef _CEKR_H
#define _CEKR_H

#include "Globals.h"
#include "SWSH.h"
#include "GKG.h"
#include "CEKG.h"
#include "FT.h"

//! Circular Equatorial Kerr Radiation Class
/*! CEKR describes methods for computing various quantities to radiation from circular, equatorial Kerr geodesics. It relies on input from instances of the classes SWSH, FT, and CEKG.  */
class CEKR {
public:
  CEKR(SWSH *swsh_in, FT *ft_in, CEKG *cekg_in);

  Complex ZI, ZH;				//!< The complex Teukolsky amplitudes

private:
  Complex rho_func();
  Complex Cnn();
  Complex Cnmbar();
  Complex Cmbarmbar();
  Complex Ann0();
  Complex Anmbar0();
  Complex Ambarmbar0();
  Complex Anmbar1();
  Complex Ambarmbar1();
  Complex Ambarmbar2();

  // Other stuff we need in lots of places.
  int l, m;
  SWSH *swsh;
  FT *ft;
  GKG gkg;
  CEKG *cekg;
  Real r, a, E, Lz;
  Real wm, pm;
  Real Delta, dr_Delta, ddr_Delta;
  Real SN_Kay, dr_SN_Kay;
  Real r_plus;

  // The above methods are used to put together ...
  Complex Zed(const Complex R,
	      const Complex dr_R,
	      const Complex ddr_R);
};

#endif

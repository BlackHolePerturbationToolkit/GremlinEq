//---------------------------------------------------------------------------
//
// $Id: CEKG.h,v 1.3 2018/08/01 22:50:14 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial Kerr geodesics.
// Scott Hughes, 11 January 1999.
//
#ifndef _CEKG_H
#define _CEKG_H

#define PROGRADE 1
#define RETROGRADE 0

#include <cmath>
#include "Globals.h"

class CEKG {
public:
  CEKG(const int orbitsense, const Real rad, const Real spin);

  Real r, a;
  Real E, Lz, Om_phi;
};
#endif

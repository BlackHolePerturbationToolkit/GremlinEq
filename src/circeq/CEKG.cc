//---------------------------------------------------------------------------
//
// $Id: CEKG.cc,v 1.4 2018/08/02 17:17:29 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial Kerr geodesics.
// Scott Hughes, 11 January 1999.
//
#include <cmath>
#include "Globals.h"
#include "CEKG.h"

CEKG::CEKG(const int orbitsense, const Real rad, const Real spin) :
  r(rad), a(spin)
{
  if (orbitsense == PROGRADE) {
    E = Kerr::Eeqpro(r, a);
    Lz = Kerr::Lzeqpro(r, a);
    Om_phi = Kerr::Omega_phi_eqpro(r, a);
  } else {
    E = Kerr::Eeqret(r, a);
    Lz = Kerr::Lzeqret(r, a);
    Om_phi = Kerr::Omega_phi_eqret(r, a);
  }
}

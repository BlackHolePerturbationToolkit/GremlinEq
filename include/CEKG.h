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

//! Circular Equatorial Kerr Geodesic Class
/*! Defines methods for computing various quantities related to circular, equatorial Kerr geodesics. */
class CEKG {
public:
  //! A constructor for the CEKG Class.
	/*!
      Initialized with orientation (prograde/retrograde), orbital radius, and spin. Circular equatorial Kerr geodesics only.
	*/
  CEKG(const int orbitsense, const Real rad, const Real spin);

  Real r; /*!< orbital radius */
  Real a; /*!< spin paramter*/
  Real E; /*!< energy */
  Real Lz; /*!< z component of angular momentum */
  Real Om_phi; /*!< Axial Frequency */
};
#endif

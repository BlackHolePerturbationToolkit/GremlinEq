//---------------------------------------------------------------------------
//
// $Id: CEDR.h,v 1.2 2018/08/04 15:55:13 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// HDF5 format data reader for circular equatorial case
//
#ifndef _CEDR_H
#define _CEDR_H

#include <hdf5_hl.h>
#include "Globals.h"

class CEDR {
public:
  CEDR(char inname[]);
  ~CEDR();
  //
  hid_t infile;
  //
  // Once the class is constructed, this will read the data set
  // for the mode pair (l,m).  Returns 1 if the file contains
  // this dataset, returns 0 if it does not.
  //
  int ReadData(const int l, const int m);
  //
  //
  // Data sets in circular equatorial files:
  //   * Properties of the orbit
  //   * Index range for the modes that are stored
  //   * SWSH eigenvalue lambda, mode amplitudes, value of the
  //     inhomogeneous solution at the orbit, mode fluxes.
  //
  Real r, a, E, Lz, Om_phi;
  //
  int lmin, lmax, mmin, mmax;
  //
  Real lambda;
  Complex ZI, ZH;
  Complex Rin, dRin, Rup, dRup;
  Real EdotI, EdotH, LzdotI, LzdotH, rdotI, rdotH;
  
private:
	bool fileExists(const std::string& filename);		
  
};
#endif

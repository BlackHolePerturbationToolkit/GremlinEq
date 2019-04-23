//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_lmode_v.cc,v 1.9 2019/04/23 16:28:19 sahughes Exp sahughes $
//
//---------------------------------------------------------------------------
//
// This code makes a table of fluxes as a function of (radius, frequency)
// for the purpose of comparing/calibrating EOB.
//
// This version: Very simple, makes a table of total flux versus r.
// Includes both energy and angular momentum flux.
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CEID.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 4) {
    cerr << "Another version of Circ_Eq_Fluxes." << endl;
    cerr << "  This version: Assumes data is evenly spaced in v." << endl;
    cerr << "  Output consists of energy for a given l mode," << endl;
    cerr << "  summed over all m for that l."
	 << endl << endl;
    cerr << "Arguments:  1. data file basename  2. Nhi  3. l" << endl;
    exit(1);
  }
  char basename[150], fluxname[150];
  FILE *fluxfile;
  //
  sprintf(basename, "%s", argv[1]);
  const int Nhi = atoi(argv[2]);
  const int l = atoi(argv[3]);
  CEID ceid(basename, Nhi, 0);
  //
  // Start generating the table.
  //
  sprintf(fluxname, "%s.flux_l%d", basename, l);
  if (!(fluxfile = fopen(fluxname, "w"))) {
    cerr << "Error opening " << fluxname << endl;
    exit(3);
  }
  //
  // Test pro/retrograde
  //
  Real rad, v, EdotHl, EdotIl, Edotnorm;
  int j;
  rad = ceid.r_arr[1];
  ceid.Get_rdot_omega(rad);
  for (j = 0; j < ceid.jmax; j++) {
    rad = ceid.r_arr[j];
    ceid.Get_rdot_omega(rad);
    v = pow(fabs(ceid.Omega), 1./3.);
    Edotnorm = (32./5.)*pow(v, 10.);
    ceid.Get_fluxes(rad);
    ceid.Get_flux_lmode(rad, l, EdotHl, EdotIl);
    //
    fprintf(fluxfile, "%lf %lf %lf %17.16e %17.16e %17.16e %17.16e %e %e\n",
	    rad, ceid.Omega, v,
	    EdotHl/Edotnorm, EdotIl/Edotnorm,
	    ceid.EdotH/Edotnorm, ceid.EdotI/Edotnorm,
	    EdotHl/ceid.EdotH, EdotIl/ceid.EdotI);
  }
  fclose(fluxfile);
}

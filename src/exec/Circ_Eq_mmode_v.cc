//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_mmode_v.cc,v 1.8 2018/08/04 16:18:40 sahughes Exp $
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
    cerr << "Arguments:  1. data file basename  2. Nhi  3. m" << endl;
    exit(1);
  }
  char basename[150], fluxname[150];
  FILE *fluxfile;
  //
  sprintf(basename, "%s", argv[1]);
  const int Nhi = atoi(argv[2]);
  const int m = atoi(argv[3]);
  CEID ceid(basename, Nhi, 0);
  //
  // Start generating the table.
  //
  sprintf(fluxname, "%s.flux_m%d", basename, m);
  if (!(fluxfile = fopen(fluxname, "w"))) {
    cerr << "Error opening " << fluxname << endl;
    exit(3);
  }
  int j;
  Real rad, v, EdotHm, EdotIm, Edotnorm;
  for (j = 1; j <= ceid.jmax; j++) {
    rad = ceid.r_arr[j];
    ceid.Get_rdot_omega(rad);
    v = pow(fabs(ceid.Omega), 1./3.);
    Edotnorm = (32./5.)*pow(v, 10.);
    ceid.Get_fluxes(rad);
    ceid.Get_flux_mmode(rad, m, EdotHm, EdotIm);
    //
    fprintf(fluxfile, "%lf %lf %lf %17.16e %17.16e %17.16e %17.16e\n",
	    rad, ceid.Omega, v,
	    EdotHm/Edotnorm, EdotIm/Edotnorm,
	    ceid.EdotH/Edotnorm, ceid.EdotI/Edotnorm);
  }
  fclose(fluxfile);
}

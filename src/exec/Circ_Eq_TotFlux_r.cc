//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_TotFlux_r.cc,v 1.10 2018/08/04 16:18:40 sahughes Exp $
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
  if(argc != 5) {
    cerr << "Original version of Circ_Eq_Fluxes." << endl;
    cerr << "  This version: Assumes data is evenly spaced in r." << endl;
    cerr << "  Output consists of summed energy and angular momentum fluxes."
	 << endl << endl;
    cerr << "Arguments:  1. data file basename" << endl;
    cerr << "            2. rmin  3. rmax  4. dr" << endl;
    exit(1);
  }
  char basename[40], fluxname[40];
  FILE *fluxfile;
  //
  sprintf(basename, "%s", argv[1]);
  const Real rmin = (Real)atof(argv[2]);
  const Real rmax = (Real)atof(argv[3]);
  const Real dr = (Real)atof(argv[4]);
  //
  CEID ceid(basename, rmin, rmax, dr, 0);
  //
  // Start generating the table.
  //
  sprintf(fluxname, "%s.flux", basename);
  if (!(fluxfile = fopen(fluxname, "w"))) {
    cerr << "Error opening " << fluxname << endl;
    exit(3);
  }
  int j;
  for (j = 1; j <= ceid.jmax; j++) {
    Real rad = ceid.r_arr[j];
    ceid.Get_rdot_omega(rad);
    ceid.Get_fluxes(rad);
    //
    fprintf(fluxfile, "%lf %lf %lf %e %e %e %e\n",
	    rad, ceid.Omega, pow(ceid.Omega, 1./3.),
	    ceid.EdotH, ceid.EdotI,
	    ceid.LzdotH, ceid.LzdotI);
  }
  fclose(fluxfile);
}

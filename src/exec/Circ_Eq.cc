//---------------------------------------------------------------------------
//
// $Id: Circ_Eq.cc,v 1.13 2018/08/03 00:22:35 sahughes Exp $
//
//---------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CETD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  //
  if(argc != 7) {
    cerr << "Arguments: 1. r  2. a" << endl;
    cerr << "           3. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           4. l  5. m  6. output file basename" << endl;
    exit(1);
  }

  const Real r = (Real)atof(argv[1]);
  const Real a = (Real)atof(argv[2]);
  const int proret = atoi(argv[3]);
  const int l = atoi(argv[4]);
  const int m = atoi(argv[5]);

  CETD cetd(proret, r, a, argv[6]);
  cetd.DoHarmonic(l, m);

  fprintf(stdout, "%13.12e %17.16e %17.16e %17.16e %17.16e\n", r,
	  cetd.ZI.real(), cetd.ZI.imag(), cetd.ZH.real(), cetd.ZH.imag());
}

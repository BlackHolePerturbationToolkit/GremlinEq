//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Seq.cc,v 1.21 2018/08/02 16:28:35 sahughes Exp $
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

  char basename[50];
  if(argc != 8) {
    cerr << "Arguments: 1. r_start  2. r_stop  3. a" << endl;
    cerr << "           4. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           5. lmax  6. dr" << endl;
    cerr << "           7. output file basename" << endl << endl;
    cerr << " NOTE: r_start assumed to be larger than r_stop." << endl;
    exit(1);
  }
  const Real r_start = (Real)atof(argv[1]);
  const Real r_stop = (Real)atof(argv[2]);
  const Real a = (Real)atof(argv[3]);
  const int proret = atoi(argv[4]);
  // Real r_isco;
  // if (proret)
  //   r_isco = Kerr::isco_pro(a);
  // else
  //   r_isco = Kerr::isco_ret(a);
  //  const Real r_min = Max(r_stop, r_isco);
  const Real r_min = r_stop;
  const int lmax = atoi(argv[5]);
  const Real dr = (Real)atof(argv[6]);
  //
  Real r;
  for (r = r_start; r >= r_min - 0.01*dr; r -= dr) {
    if (dr < 0.01)
      sprintf(basename, "%s_r%4.3lf", argv[7], r);
    else if (dr < 0.1)
      sprintf(basename, "%s_r%3.2lf", argv[7], r);
    else
      sprintf(basename, "%s_r%2.1lf", argv[7], r);
    //
    CETD cetd(proret, r, a, basename);
    //
    cetd.Driver(lmax);
  }
}

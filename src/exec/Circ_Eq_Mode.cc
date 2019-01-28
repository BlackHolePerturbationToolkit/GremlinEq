//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Mode.cc,v 1.6 2019/01/20 17:54:16 sahughes Exp $
//
//---------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "Globals.h"
#include "CEKG.h"
#include "SWSH.h"
#include "FT.h"
#include "CEKR.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  //
  if(argc != 6) {
    cerr << "Arguments: 1. a  2. r  3. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           4. l  5. m" << endl;
    exit(1);
  }
  const Real a = (Real)atof(argv[1]);
  const Real r = (Real)atof(argv[2]);
  const int proret = atoi(argv[3]);
  const int l = atoi(argv[4]);
  const int m = atoi(argv[5]);

  CEKG cekg(m, proret, r, a);
  SWSH swsh(-2, l, m, a*cekg.wm);
  FT ft(l, m, r, a, cekg.wm, swsh.lambda, 1.e-12);
  CEKR cekr(&swsh, &ft, &cekg);

  fprintf(stdout, "Re(ZH) = %e\n", cekr.ZH.real());
  fprintf(stdout, "Im(ZH) = %e\n", cekr.ZH.imag());
  fprintf(stdout, "Re(ZI) = %e\n", cekr.ZI.real());
  fprintf(stdout, "Im(ZI) = %e\n", cekr.ZI.imag());
}

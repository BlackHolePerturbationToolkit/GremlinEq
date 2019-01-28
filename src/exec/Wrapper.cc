//---------------------------------------------------------------------------
//
// $Id: Wrapper.cc,v 1.18 2019/01/27 09:57:28 sahughes Exp $
//
//---------------------------------------------------------------------------
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CEKG.h"
#include "SWSH.h"
#include "FT.h"
#include "CEKR.h"
#include "CETD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  if(argc != 5) {
    cerr << "Arguments: 1. l  2. m  3. a  4. r" << endl;
    exit(1);
  }
  const int l = atoi(argv[1]);
  const int m = atoi(argv[2]);
  const Real a = (Real)atof(argv[3]);
  const Real r = (Real)atof(argv[4]);

  CEKG cekg(PROGRADE, r, a);
  SWSH swsh(-2, l, m, a*m*cekg.Om_phi);
  Real x;

  cout << a*m*cekg.Om_phi << " " << swsh.N << endl;

//   for (x = -0.99999; x <= 1.00000001; x += 0.0001) {
//     const Real S = swsh.spheroid(x);
//     const Real l2S = swsh.l2dagspheroid(x);
//     const Real l1l2S = swsh.l1dagl2dagspheroid(x);
//     fprintf(stdout, "%lf %e %e %e\n", x, S, l2S, l1l2S);
//   }
}


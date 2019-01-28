//---------------------------------------------------------------------------
//
// $Id: TidalHEq.cc,v 1.17 2019/01/20 17:54:35 sahughes Exp $
//
//---------------------------------------------------------------------------

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
#include "TidalH.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  if(argc != 4 && argc != 5) {
    cerr << "Arguments: 1. a  2. r  3. l value to focus on" << endl;
    cerr << "           4. m value to focus on" << endl;
    cerr << " Argument 4 is optional; if it is not included, the sum" << endl;
    cerr << " is done over -l <= m < = l.  Otherwise, it is only done" << endl;
    cerr << " for the m that you provide.  (Actually, adding m and -m" << endl;
    cerr << " for m != 0)." << endl << endl;
    exit(1);
  }
  const Real a = (Real)atof(argv[1]);
  const Real r = (Real)atof(argv[2]);
  const int lp = atoi(argv[3]);
  int m, mhere;
  if (argc == 5)
    mhere = abs(atoi(argv[4]));

  const Real photorb = 2.*(1. + cos(2.*acos(-a)/3.));

  if (fabs(a) > 1.) {
    cerr << "Spin greater than 1" << endl;
    exit(0);
  }
  if (r < photorb) {
    cerr << "Radius is inside photon orbit!" << endl;
    exit(0);
  }
  if (r < Kerr::isco_pro(a))
    cerr << "Warning: radius inside the ISCO." << endl;

  TidalH th(a, lp);

  Real x, psi, R1;
  int l;
  for (l = 2; l <= lp; l ++) {
    for (m = -l; m <= l; m++) {
      CEKG cekg(1, r, a);
      SWSH swshneg(-2, l, m, a*m*cekg.Om_phi);
      FT ft(l, m, r, a, m*cekg.Om_phi, swshneg.lambda, 3.e-14);
      CEKR cekr(&swshneg, &ft, &cekg);
      //
      SWSH swshpos(2, l, m, a*m*cekg.Om_phi);
      //
      // Load data in th.
      //
      th.swshp[l][m] = swshpos;
      th.w[m] = m*cekg.Om_phi;
      th.p[m] = m*cekg.Om_phi - m*a/(2.*Kerr::rplus(a));
      th.lamb[l][m] = swshneg.lambda;
      th.ZH[l][m] = cekr.ZH*r*r*r;
    }
  }
  th.loadClm();
  x = 0.;
  for (psi = 0.0; psi <= 2.0000001*M_PI; psi += 0.00001*M_PI) {
    R1 = 0.;
    for (m = -lp; m <= lp; m++) {
      if (argc != 5 || abs(m) == mhere)
	R1 += th.R1lm(lp, m, x, psi);
    }
    fprintf(stdout, "%e %e\n", psi, R1);
  }
}

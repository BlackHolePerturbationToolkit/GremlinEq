//---------------------------------------------------------------------------
//
// $Id: EmbedEq.cc,v 1.18 2019/01/20 17:54:16 sahughes Exp $
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
#include "TidalH.h"

#define JMAX 12
#define JMAXP (JMAX+1)
#define K 5
#define TINY 1.0e-25
#define EPS 2.e-15

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  if(argc != 4) {
    cerr << "Arguments: 1. a  2. r  3. lmax" << endl;
    exit(1);
  }
  const Real a = (Real)atof(argv[1]);
  const Real r = (Real)atof(argv[2]);
  const int lmax = atoi(argv[3]);

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

  const int lmin = 2;
  int l, m;
  TidalH th(a, lmax);
  //
  // Compute the angle at which the orbiting body sits in the co-rotating
  // frame.
  //
  Real rs, rb;
  const Real rp = Kerr::rplus(a);
  const Real rm = Kerr::rminus(a);
  const Real Omorb = Kerr::Omega_phi_eqpro(r, a);
  const Real soma2 = sqrt((1. - a)*(1. + a));
  if (fabs(a) > 1.e-8) {
    rs = r + (rp/soma2)*log(r/rp - 1.) - (rm/soma2)*log(r/rm - 1.);
    rb = 0.5*(a/soma2)*log((r - rp)/(r - rm));
  } else {
    rs = r + 2.*log(0.5*r - 1.);
    rb = 0.5*a*log(1. - 2./r);
  }
  const Real Deltapsi_orb = rb - Omorb*rs;
  for (l = lmin; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      CEKG cekg(1, r, a);
      SWSH swshneg(-2, l, m, a*m*cekg.Om_phi);
      FT ft(l, m, r, a, m*cekg.Om_phi, swshneg.lambda, 3.e-14);
      CEKR cekr(&swshneg, &ft, &cekg);
      //
      SWSH swshpos(2, l, m, a*m*cekg.Om_phi);
      //
      // Load data into th.
      //
      th.swshp[l][m] = swshpos;
      th.w[m] = m*cekg.Om_phi;
      th.p[m] = m*cekg.Om_phi - m*a/(2.*rp);
      th.lamb[l][m] = swshneg.lambda;
      th.ZH[l][m] = cekr.ZH*r*r*r; // Notice r^3 scaled out.
    }
  }
  th.loadClm();
  th.loadEpslm();
  //
  // Now produce some output.
  //
  Real x = 0, psi;
  Real R1;
  Real radmax = 0.;
  Real radorb;
  Real psi_of_radmax;
  Real epsilon_rad;
  int ell, ellmin = lmin, ellmax = lmax;
  for (psi = Deltapsi_orb; psi <= Deltapsi_orb + 2.0000001*M_PI;
       psi += 0.0001*M_PI) {
    R1 = 0.;
    epsilon_rad = 0.;
    for (ell = ellmin; ell <= ellmax; ell++) {
      for (m = -ell; m <= ell; m++) {
	R1 += th.R1lm(ell, m, x, psi);
  	epsilon_rad += th.epsr(ell, m, x, psi);
      }
    }
    Real scalefactor = 15.;
    Real rad0 = 2.;
    Real rad = rad0 + epsilon_rad/scalefactor;

    cout << rad0*cos(psi) << " " << rad0*sin(psi) << " "
  	 << rad*cos(psi) << " " << rad*sin(psi) << " "
  	 << psi << " " << R1 << " "
	 << epsilon_rad
	 << endl;

    if (psi < 0.5*M_PI || psi > 1.5*M_PI) {
      if (rad > radmax) {
	radmax = rad;
	psi_of_radmax = psi;
      }
    }
    if (fabs(psi - Deltapsi_orb) < 1.e-10)
      radorb = rad;
  }

  if (psi_of_radmax > M_PI) psi_of_radmax -= 2.*M_PI;

  cout << endl << "# Orbit position: " << radorb
       << " at psi = " << Deltapsi_orb
       << " = " << Deltapsi_orb*180./M_PI << " degrees" << endl;
  cout << "# Coordinates: x = " << radorb*cos(Deltapsi_orb)
       << ", y = " << radorb*sin(Deltapsi_orb) << endl;
  cout << endl << "# Largest radius: " << radmax
       << " at psi = " << psi_of_radmax
       << " = " << psi_of_radmax*180./M_PI << " degrees" << endl;
  cout << "# Coordinates: x = " << radmax*cos(psi_of_radmax)
       << ", y = " << radmax*sin(psi_of_radmax) << endl;
  cout << endl << "# Opening angle of wedge: "
       << (psi_of_radmax - Deltapsi_orb)*180./M_PI << " degrees" << endl;
}

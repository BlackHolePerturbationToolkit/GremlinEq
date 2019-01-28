//---------------------------------------------------------------------------
//
// $Id: EmbedSurf.cc,v 1.12 2019/01/27 00:50:57 sahughes Exp $
//
//---------------------------------------------------------------------------
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_integration.h>
#include "Globals.h"
#include "CEKG.h"
#include "SWSH.h"
#include "FT.h"
#include "CEKR.h"
#include "TidalH.h"
#include "Tensors.h"

Real Gigrand(Real x, void *params)
{
  const Real a = *(Real *)params;
  const Real x2 = x*x;
  const Real x4 = x2*x2;
  const Real rp = Kerr::rplus(a);
  const Real rp2 = rp*rp;
  const Real rp4 = rp2*rp2;
  const Real rp2pa2x2 = rp2 + a*a*x2;

  const Real Hsqr = (rp4*rp4) - 6.*a*a*a*a*rp4*x2 -
    4.*a*a*a*a*a*a*rp2*x2*(1. + x2) - a*a*a*a*a*a*a*a*x2*(1. + x2 + x4);

  const Real anssqr = Hsqr/(rp2pa2x2*rp2pa2x2*rp2pa2x2);

  return(sqrt(anssqr));
}

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  if(argc != 4) {
    cerr << "Arguments: 1. a  2. r  3. lmax" << endl;
    exit(1);
  }
  Real a = (Real)atof(argv[1]);
  if (a > 0.86602540378444) {
    cerr << "This code can't do embeddings for a > sqrt(3)/2" << endl;
    exit(2);
  }
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
      th.p[m] = m*cekg.Om_phi - m*a/(2.*Kerr::rplus(a));
      th.lamb[l][m] = swshneg.lambda;
      th.ZH[l][m] = cekr.ZH*r*r*r; // Notice r^3 scaled out.
    }
  }
  th.loadClm();
  th.loadEpslm();
  //
  // Now produce some output.
  //
  Real x, psi;
  Real R1;
  Real radmax = 0.;
  Real epsilon_rad;
  Real rp = Kerr::rplus(a);
  //
  // eta and beta are defined in Smarr 1975.
  //
  Real eta = sqrt(2.*rp);
  Real beta = a/eta;
  int ell, ellmin = lmin, ellmax = lmax;
  for (x = -1.; x <= 1.000000001; x += 0.05) {
    Real xh;
    if (x <= -1.0) {
      x = -1.0;
      xh = -0.9999999;
    } else if (x >= 1.0) {
      x = 1.0;
      xh = 0.9999999;
    } else
      xh = x;
    const Real f = (1. - xh)*(1. + xh)/((1. - beta*beta*(1. - xh)*(1. + xh)));
    const Real rho0 = eta*sqrt(f);
    //
    gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(18);
    size_t nevals;
    gsl_function Igrand;
    Real Z0;
    Igrand.function = &Gigrand;
    Igrand.params = &a;
    gsl_integration_romberg(&Igrand, 0, xh, 0, 1e-14, &Z0, &nevals, w);
    //
    Real rad0 = sqrt(rho0*rho0 + Z0*Z0);
    gsl_integration_romberg_free(w);                                                                
    for (psi = 0.0; psi <= 2.0000001*M_PI; psi += 0.05*M_PI) {
      Real X0 = rho0*cos(psi);
      Real Y0 = rho0*sin(psi);
      R1 = 0.;
      epsilon_rad = 0.;
      Real scalefactor = 30.;
      for (ell = ellmin; ell <= ellmax; ell++) {
	for (m = -ell; m <= ell; m++) {
	  R1 += th.R1lm(ell, m, x, psi);
	  epsilon_rad += th.epsr(ell, m, xh, psi);
	}
      }
      Real dr = rp*epsilon_rad/scalefactor;
      
      Real dX = dr*sqrt((1. - x)*(1. + x))*cos(psi);
      Real dY = dr*sqrt((1. - x)*(1. + x))*sin(psi);
      Real dZ = dr*x;
      
      cout << X0 << " " << Y0 << " " << Z0 << " "
	   << X0 + dX << " " << Y0 + dY << " " << Z0 + dZ << " "
	   << endl;
    }
    cout << endl;
  }
}


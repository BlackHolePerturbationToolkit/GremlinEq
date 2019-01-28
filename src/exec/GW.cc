//---------------------------------------------------------------------------
//
// $Id: GW.cc,v 1.9 2019/01/20 17:54:35 sahughes Exp $
//
//---------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CEDR.h"
#include "SWSH.h"
#include "RRGW.h"
#include "Tensors.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  int i, l, m;
  const int lmax = 30;
  Real r, a, t, dhp, dhc, hp_geom, hc_geom;
  Real *w;
  //
  Complex **ZI, **ZH;
  Complex psi4, dpsi4;
  Real **Spheroid;
  RRGW rrgw;
  char inname[120], wavename[120];
  FILE *wavefile;

  if(argc != 6) {
    cerr << "Arguments: 1. Input file basename.  2. cos(theta)." << endl;
    cerr << "           3. phi (degrees)  4. dt  5. Nsteps" << endl;
    exit(1);
  }
  sprintf(inname, "%s.h5", argv[1]);
  sprintf(wavename, "%s.wave", argv[1]);
  const Real ct = (Real)atof(argv[2]);
  const Real phi = M_PI*(Real)atof(argv[3])/180.;
  const Real dt = (Real)atof(argv[4]);
  const int Nsteps = atoi(argv[5]);
  //
  // Allocate memory.
  //
  w = Tensor<Real>::vector(-lmax, lmax);
  ZI = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);
  ZH = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);
  Spheroid = Tensor<Real>::matrix(2, lmax, -lmax, lmax);
  //
  // Pre-load these guys with zeros; this way there is no
  // blank memory for the multipoles that don't get computed.
  //
  for(l = 2; l <= lmax; l++) {
    for(m = -l; m <= l; m++) {
      ZI[l][m] = Complex(0., 0.);
      ZH[l][m] = Complex(0., 0.);
    }
  }
  //
  // Open data.
  //
  CEDR cedr(inname);
  r = cedr.r; a = cedr.a;
  for (m = -lmax; m <= lmax; m++) {
    w[m] = m*cedr.Om_phi;
  }
  for (l = cedr.lmin; l <= cedr.lmax; l++) {
    for (m = cedr.mmin; m <= cedr.mmax; m++) {
      if (cedr.ReadData(l, m)) {
	ZI[l][m] = cedr.ZI; ZH[l][m] = cedr.ZH;
	Real sign = (l%2) ? -1. : 1.;
	ZH[l][-m] = sign*conj(ZH[l][m]);
	ZI[l][-m] = sign*conj(ZI[l][m]);
      }
    }
  }
  //
  // Write the waveform.
  //
  // First, do all the spheroidal harmonics.  They need only be
  // calculated once and stored.
  //
  for(l = cedr.lmin; l <= cedr.lmax; l++) {
    for(m = -l; m <= l; m++) {
      if(w[m] != 0.0) {
	SWSH swsh(-2, l, m, a*w[m]);
	Spheroid[l][m] = swsh.spheroid(ct);
      }
      else 
	Spheroid[l][m] = 0.0;
    }
  }
  //
  wavefile = fopen(wavename, "w");
  for(i = 0; i <= Nsteps; i++) {
    t = ((Real)i)*dt;
    hp_geom = 0.; hc_geom = 0.;
    psi4 = Complex(0., 0.);
    for(l = cedr.lmin; l <= cedr.lmax; l++) {
      for(m = -l; m <= l; m++) {
	if(w[m] != 0) {
	  rrgw.Wave(m, t, phi, Spheroid[l][m], w[m], ZI[l][m], dhp, dhc);
	  rrgw.Psi4(m, t, phi, Spheroid[l][m], w[m], ZI[l][m], dpsi4);
	  hp_geom += dhp; hc_geom += dhc;
	  psi4 += dpsi4;
	}
      }
    }
    fprintf(wavefile, "%e %e %e %e %e %e\n", t, hp_geom, hc_geom,
	    psi4.real(), psi4.imag(), abs(psi4));
  }
  fclose(wavefile);
  //
  // Free memory.
  //
  Tensor<Real>::free_vector(w, -lmax, lmax);
  Tensor<Complex>::free_matrix(ZI, 2, lmax, -lmax, lmax);
  Tensor<Complex>::free_matrix(ZH, 2, lmax, -lmax, lmax);
  Tensor<Real>::free_matrix(Spheroid, 2, lmax, -lmax, lmax);
}

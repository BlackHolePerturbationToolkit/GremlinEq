//---------------------------------------------------------------------------
//
// $Id: Fluxes.cc,v 1.12 2018/08/04 16:19:14 sahughes Exp $
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
#include "Tensors.h"
#include "SWSH.h"
#include "RRGW.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  int i, l, m, harm_count = 0;
  const int lrange = 100;
  Real E, Lz, Q;
  Real r, a, Om_phi;
  //
  // Word of caution: Edot and Lzdot in the datafile are the
  // energy and angular momentum rate of change in the radiation.
  // The change in the particle's energy and angular momentum is
  // Edot_particle = -Edot, Lzdot_particle = -Lzdot.
  //
  Complex **ZI, **ZH;
  Real **EdotI, **EdotH;
  Real **LzdotI, **LzdotH;
  Real **rdotI, **rdotH;
  Real EdotI_tot = 0., LzdotI_tot = 0., rdotI_tot = 0.;
  Real EdotH_tot = 0., LzdotH_tot = 0., rdotH_tot = 0.;
  Real Edot_tot = 0., Lzdot_tot = 0., rdot = 0.;
  RRGW rrgw;
  char inname[100], modeset[50];

  if(argc != 2) {
    cerr << "Arguments: 1. Input filename." << endl;
    exit(1);
  }
  sprintf(inname, "%s", argv[1]);
  //
  // Allocate memory.
  //
  ZI = Tensor<Complex>::matrix(2, lrange, -lrange, lrange);
  ZH = Tensor<Complex>::matrix(2, lrange, -lrange, lrange);
  EdotI = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  LzdotI = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  rdotI = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  EdotH = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  LzdotH = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  rdotH = Tensor<Real>::matrix(2, lrange, -lrange, lrange);
  //
  // Pre-load these guys with zeros; this way there is no
  // blank memory for the multipoles that don't get computed.
  //
  for(l = 2; l <= lrange; l++) {
    for(m = -l; m <= l; m++) {
      ZI[l][m] = Complex(0., 0.); ZH[l][m] = Complex(0., 0.);
      EdotI[l][m] = 0.; LzdotI[l][m] = 0.; rdotI[l][m] = 0.;
      EdotH[l][m] = 0.; LzdotH[l][m] = 0.; rdotH[l][m] = 0.;
    }
  }
  //
  // Open data
  //
  CEDR cedr(inname);
  if(!cedr.exists){
	  cout << "File " << inname << " does not exist" << endl;
	  exit(0);
  }
  
  r = cedr.r; a = cedr.a;
  E = cedr.E; Lz = cedr.Lz;
  Om_phi = cedr.Om_phi;
  for (l = cedr.lmin; l <= cedr.lmax; l++) {
    for (m = cedr.mmin; m <= cedr.mmax; m++) {
      if (cedr.ReadData(l, m)) {
	ZI[l][m] = cedr.ZI; ZH[l][m] = cedr.ZH;
	EdotI[l][m] = cedr.EdotI; EdotH[l][m] = cedr.EdotH;
	LzdotI[l][m] = cedr.LzdotI; LzdotH[l][m] = cedr.LzdotH;
	rdotI[l][m] = cedr.rdotI; rdotH[l][m] = cedr.rdotH;
	harm_count++;
      }
    }
  }
  fprintf(stderr,"%d harmonics computed for this orbit.\n", harm_count);
  //
  // Use symmetry to fill in blank spots.
  //
  for(l = 2; l <= lrange; l++) {
    for(m = -l; m <= l; m++) {
      const Real sign = (l%2 ? -1 : 1.);
      if (ZI[l][m].real() == 0. && ZI[l][m].imag() == 0.)
	ZI[l][m] = sign*conj(ZI[l][-m]);
      if (ZH[l][m].real() == 0. && ZH[l][m].imag() == 0.)
	ZH[l][m] = sign*conj(ZH[l][-m]);
      if (EdotI[l][m] == 0.) EdotI[l][m] = EdotI[l][-m];
      if (LzdotI[l][m] == 0.) LzdotI[l][m] = LzdotI[l][-m];
      if (rdotI[l][m] == 0.) rdotI[l][m] = rdotI[l][-m];
      if (EdotH[l][m] == 0.) EdotH[l][m] = EdotH[l][-m];
      if (LzdotH[l][m] == 0.) LzdotH[l][m] = LzdotH[l][-m];
      if (rdotH[l][m] == 0.) rdotH[l][m] = rdotH[l][-m];
    }
  }
  for(l = 2; l <= lrange; l++) {
    for(m = -l; m <= l; m++) {
      EdotI_tot += EdotI[l][m];
      LzdotI_tot += LzdotI[l][m];
      rdotI_tot += rdotI[l][m];
      EdotH_tot += EdotH[l][m];
      LzdotH_tot += LzdotH[l][m];
      rdotH_tot += rdotH[l][m];
    }
  }
  Edot_tot = EdotI_tot + EdotH_tot;
  Lzdot_tot = LzdotI_tot + LzdotH_tot;
  rdot = rdotI_tot + rdotH_tot;
  //
  Real v = pow(fabs(Om_phi), 1./3.);
  Real EdotN = (32./5.)*pow(v, 10.);

  fprintf(stdout, "EdotH: %.16e EdotI: %.16e Edottot: %.16e\n",
 	  EdotH_tot/EdotN, EdotI_tot/EdotN, (EdotH_tot + EdotI_tot)/EdotN);
  //
  // Output all the harmonics
  //
  for (l = cedr.lmin; l <= cedr.lmax; l++) {
    for (m = -l; m <= l; m++) {
      if (EdotI[l][m] + EdotH[l][m] != 0.0)
	fprintf(stdout,
		"%d %d %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
		l, m,
		ZI[l][m].real(), ZI[l][m].imag(),
		ZH[l][m].real(), ZH[l][m].imag(),
		EdotH[l][m], EdotI[l][m],
		LzdotH[l][m], LzdotI[l][m]);
    }
  }
  //
  // Free memory.
  //
  Tensor<Complex>::free_matrix(ZI, 2, lrange, -lrange, lrange);
  Tensor<Complex>::free_matrix(ZH, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(EdotI, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(LzdotI, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(rdotI, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(EdotH, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(LzdotH, 2, lrange, -lrange, lrange);
  Tensor<Real>::free_matrix(rdotH, 2, lrange, -lrange, lrange);
}

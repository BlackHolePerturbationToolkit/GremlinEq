//---------------------------------------------------------------------------
//
// $Id: TidalHSurf.cc,v 1.12 2019/01/20 17:54:35 sahughes Exp $
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
#include "Tensors.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  if(argc != 3) {
    cerr << "Arguments: 1. r  2. a" << endl;
    exit(1);
  }
  const Real r = (Real)atof(argv[1]);
  const Real a = (Real)atof(argv[2]);

  // Globally true for everything
  const Real rp = Kerr::rplus(a);
  const Real wH = a/(2.*rp);
  const Real eps = sqrt((1. - a)*(1. + a))/(4.*rp);

  int lmax = 12;
  int l, m;
 
  Real x, phi;
  Real R1;

  Complex **ZH;
  ZH = Tensor<Complex>::matrix(2, lmax, -lmax, lmax);

  Real *w; w = Tensor<Real>::vector(-lmax, lmax);
  Real *p; p = Tensor<Real>::vector(-lmax, lmax);
  Real **lamb; lamb = Tensor<Real>::matrix(2, lmax, -lmax, lmax);

  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      CEKG cekg(1, r, a);
      SWSH swshneg(-2, l, m, a*m*cekg.Om_phi);
      FT ft(l, m, r, a, m*cekg.Om_phi, swshneg.lambda, 3.e-14);
      CEKR cekr(&swshneg, &ft, &cekg);
      ZH[l][m] = cekr.ZH*r*r*r;
      w[m] = m*cekg.Om_phi;
      p[m] = w[m] - m*wH; 
      lamb[l][m] = swshneg.lambda;
    }
  }

  for (x = -1.; x <= 1.0000001; x += 0.02) {
    for (phi = 0.0; phi <= 2.0000001*M_PI; phi += 0.02*M_PI) {
      R1 = 0.;
      for (l = 2; l <= lmax; l++) {
	for (m = -l; m <= l; m++) {
	  SWSH swshpos(2, l, m, a*w[m]);
	  //
	  // Some auxiliary quantities we need
	  //
	  const Real lp2 = lamb[l][m] + 2.;
	  const Real tmp1 = lp2*lp2 + 4.*a*w[m]*(((Real)m) - a*w[m]);
	  const Real tmp2 = lamb[l][m]*lamb[l][m] +
	    36.*a*w[m]*(((Real)m) - a*w[m]);
	  const Real tmp3 = 48.*a*w[m]*(2.*lamb[l][m] + 3.)*
	    (2.*a*w[m] - ((Real)m));
	  
	  Real abs_sqr_clm = tmp1*tmp2 + tmp3 + 144.*w[m]*w[m]*(1. - a*a);
	  const Real Im_clm = 12.*w[m];
	  const Real Re_clm = sqrt(abs_sqr_clm - Im_clm*Im_clm);
	  const Complex clm = Re_clm + II*Im_clm;
	  
	  Real xhere;
	  if (x <= -1.0) {
	    x = -1.0;
	    xhere = -0.9999999;
	  } else if (x >= 1.0) {
	    x = 1.0;
	    xhere = 0.9999999;
	  } else
	    xhere = x;
	  
	  const Complex R1lm = 512.*rp*rp*(p[m] + 4.*II*eps)*exp(II*m*phi)*
	    (II*p[m] - 2.*eps)*ZH[l][m]*
	    swshpos.edthbaredthbarspheroid(a, xhere)/clm;
	  
	  R1 += R1lm.imag();
	}
      }
      Real scalefactor = 16.;
      Real Gauss1 = 0.5*R1/scalefactor;
      Real rad0 = 2.;
      Real rad = rad0 - 0.5*Gauss1;
      Real XX = rad*sqrt((1. - x)*(1. + x))*cos(phi);
      Real YY = rad*sqrt((1. - x)*(1. + x))*sin(phi);
      Real ZZ = rad*x;
      Real xx = rad0*sqrt((1. - x)*(1. + x))*cos(phi);
      Real yy = rad0*sqrt((1. - x)*(1. + x))*sin(phi);
      Real zz = rad0*x;

      cout << XX << " " << YY << " " << ZZ << " "
	   << xx << " " << yy << " " << zz << " "
	   << x << " " << phi << " " << R1 << endl;
    }
    cout << endl;
  }
}


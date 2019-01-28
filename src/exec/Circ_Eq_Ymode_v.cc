//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Ymode_v.cc,v 1.17 2018/08/04 16:18:40 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// This code makes a table of fluxes as a function of (radius, frequency, v)
// for the purpose of comparing/calibrating EOB.
//
// This version: Makes a table of flux versus v.  Provides energy flux
// for a particular harmonic, but rotated into spherical harmonics
// rather than spheroidal harmonics.  Also gives wave amplitude hlm.
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CEID.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 7) {
    cerr << "Another version for working with data parameterized by v."
	 << endl;
    cerr << "  This version: Assumes data evenly spaced in v," << endl;
    cerr << "  with Nhi total files." << endl;
    cerr << "  Output consists of energy flux AND waveform for a single"
	 << endl;
    cerr << "  *spherical* harmonic mode (l,m).  Data up to lmax are" << endl;
    cerr << "  included in the calculation.  Buffer is related to" << endl;
    cerr << "  the number of spheroidal ls that are used to convert" << endl;
    cerr << "  to spherical l." << endl << endl;
    cerr << "Arguments:  1. data file basename  2. Nhi" << endl;
    cerr << "            3. l  4. m  5. lmax  6. buffer" << endl;
    cerr << endl;
    cerr << " Spheroidal* data at indices up to" << endl;
    cerr << " l = lmax + buffer will be included in the" << endl;
    cerr << " conversion to the spherical harmonic basis." << endl;
    exit(1);
  }
  char basename[140], fluxname[140];
  FILE *fluxfile;
  //
  sprintf(basename, "%s", argv[1]);
  const int Nhi = atoi(argv[2]);
  const int l = atoi(argv[3]);
  const int m = atoi(argv[4]);
  const int lmax = atoi(argv[5]);
  const int BUFFER = atoi(argv[6]);
  //
  CEID ceid(basename, Nhi, lmax + BUFFER);
  //
  // Start generating the table.
  //
  sprintf(fluxname, "%s.Yflux_%d%d", basename, l, m);
  if (!(fluxfile = fopen(fluxname, "w"))) {
    cerr << "Error opening " << fluxname << endl;
    exit(3);
  }
  int i;
  Complex CIlm, hlm;
  Real EdotI_Y, EdotH_Y;
  ceid.Spline_CIlm(l, m);
  for (i = Nhi; i >= 1; i--) {
    Real rad = ceid.r_arr[i];
    ceid.Get_rdot_omega(rad);
    Real v = pow(fabs(ceid.Omega), 1./3.);
    Real Edotnorm = (32./5.)*pow(v, 10.);
    Real q = ceid.a;
    Real EPNT_Spin_Less = 1. - 3.71130952380952381*pow(v,2.)
      + 12.5663706143591730*pow(v,3.) - 4.92846119929453263*pow(v,4.)
      - 38.2928354546934470*pow(v,5.)
      + (115.731716675611331 - 16.3047619047619048*log(v))*pow(v,6.)
      - 101.509595959741638*pow(v,7.)
      - (117.504390722677329 - 52.7430839002267574*log(v))*pow(v,8.)
      + (719.128342233429702 - 204.891680874122896*log(v))*pow(v,9.)
      - (1216.90699131704196 - 116.639876594109398*log(v))*pow(v,10.)
      + (958.934970119566876 + 473.624478174230623*log(v))*pow(v,11.);
    Real EPNT_Spin_Full = -11./4.*q*pow(v,3.) + 33./16.*pow(q,2.)*pow(v,4.)
      - 59./16.*q*pow(v,5.)
      - (65./6.*M_PI*q - 611./504.*pow(q,2.))*pow(v,6.)
      + (162035./3888.*q + 65./8.*M_PI*pow(q,2.)
    	 - 71./24.*pow(q,3.))*pow(v,7.)
      + (-359./14.*M_PI*q + 22667./4536.*pow(q,2.)
    	 + 17./16.*pow(q,4.))*pow(v,8.);  
    Real EPNT = EPNT_Spin_Less + EPNT_Spin_Full;  

    ceid.Get_flux_Ymode(rad, l, m, EdotH_Y, EdotI_Y);
    ceid.Get_CIlm(rad, 0., l, m, CIlm);
    hlm = -2.*CIlm/(m*m*ceid.Omega*ceid.Omega);
    EdotH_Y *= 2.; // want both m and -m ...
    EdotI_Y *= 2.; // want both m and -m ...
    fprintf(fluxfile,
	    "%13.12e %13.12e %17.16e %17.16e %17.16e %17.16e %17.16e\n",
	    v, rad, EdotH_Y/Edotnorm, EdotI_Y/Edotnorm,
	    (EdotH_Y + EdotI_Y)/Edotnorm,
	    hlm.real(), hlm.imag());
  }
  fclose(fluxfile);
}

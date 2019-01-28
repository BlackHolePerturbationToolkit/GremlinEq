//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_TotFlux_v.cc,v 1.19 2018/08/04 16:18:40 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// This code makes a table of fluxes as a function of (radius, frequency, v)
// for the purpose of comparing/calibrating EOB.
//
// This version: Makes a table of flux versus v.  Does only energy
// flux summed over modes.  Gives flux to infinity, horizon, and
// total; also compares with PN expansion.
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
  if(argc != 3) {
    cerr << "Version for working with data parameterized by v." << endl;
    cerr << "  This version: Assumes data is evenly spaced in v,"<< endl;
    cerr << "  with Nhi total files." << endl;
    cerr << "  Output consists of summed energy flux and some diagnostics."
	 << endl << endl;
    cerr << "Arguments:  1. data file basename  2. Nhi" << endl;
    exit(1);
  }
  char basename[40], fluxname[40];
  FILE *fluxfile;
  //
  sprintf(basename, "%s", argv[1]);
  const int Nhi = atoi(argv[2]);
  //
  CEID ceid(basename, Nhi, 0);
  //
  // Start generating the table.
  //
  sprintf(fluxname, "%s.flux", basename);
  if (!(fluxfile = fopen(fluxname, "w"))) {
    cerr << "Error opening " << fluxname << endl;
    exit(3);
  }
  int i;
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

    ceid.Get_fluxes(rad);

    // fprintf(fluxfile, "%13.12e %13.12e %17.16e %17.16e %17.16e %e %d\n",
    // 	    v, rad,
    // 	    ceid.EdotH/Edotnorm, ceid.EdotI/Edotnorm,
    // 	    (ceid.EdotH + ceid.EdotI)/Edotnorm,
    // 	    EPNT - ceid.EdotI/Edotnorm,
    // 	    ceid.max_l_computed[i]);
    int prograde;
    if (ceid.Omega >= 0) prograde = 1;
    else prograde = 0;
    Real en, lz;
    if (prograde) {
      en = Kerr::Eeqpro(rad, q);
      lz = Kerr::Lzeqpro(rad, q);
    } else {
      en = Kerr::Eeqret(rad, q);
      lz = Kerr::Lzeqret(rad, q);
    }
    Real delta = Kerr::Delta(rad, q);
    Real r2pa2 = rad*rad + q*q;
    Real Gamma = en*(r2pa2*r2pa2/delta - q*q) - 2.*rad*q*lz/delta;

    // fprintf(fluxfile, "%.8e %.8e %.10e %.10e %.10e %.10e\n",
    // 	    v, rad,
    // 	    ceid.rdot,
    // 	    ceid.EdotH + ceid.EdotI,
    // 	    (ceid.EdotH + ceid.EdotI)/Edotnorm,
    // 	    (ceid.LzdotH + ceid.LzdotI)*Gamma);
    const Real mrsq = 1.e-8;
    fprintf(fluxfile, "%.8e %.8e %.8e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
  	    v, rad, ceid.Omega,
	    ceid.rdot,
	    ceid.rdot_noH,
  	    mrsq*(ceid.EdotH + ceid.EdotI),
  	    mrsq*(ceid.EdotI),
  	    (ceid.EdotH + ceid.EdotI)/Edotnorm,
  	    ceid.EdotI/Edotnorm,
  	    (ceid.LzdotH + ceid.LzdotI)*Gamma,
  	    ceid.LzdotI*Gamma);
  }
  fclose(fluxfile);
}

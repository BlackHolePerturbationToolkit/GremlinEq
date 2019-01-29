//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Clm2.cc,v 1.26 2018/08/04 16:18:40 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// This code will perform interpolations along a radiation reaction
// sequence to generate a trajectory through the parameter space of
// circular, equatorial Kerr geodesic orbits.
//
// Input parameters: just needs masses of the bodies and the filename.
// The initial radius is in the file, and the masses then set all the
// timescales.
//
// In particular, note that the time used in the integrator is related
// to Boyer-Lindquist time by t_BL = (M*M/mu) t_Code.  We actually
// output everything per unit mass, so we output \hat t_BL = (M/mu)
// t_Code = t_Code/eta.  In the final computation and output, this
// scaling is correctly accounted for in all quantities.
//
// This version makes the quantities C_lm used in numerical relativity,
// which is basically the coefficients of psi_4 on a particular spin-
// weighted spherical harmonic.  ONLY WORKS FOR SCHWARZSCHILD RIGHT NOW.
// If used for Kerr, it will give the projection on a spin-weighted
// *spheroidal* harmonic.
//
// Restriction to Schwarzschild only lifted, 8 June 2010.
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CEID.h"

//extern "C" int isnan(double);

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 12) {
    cerr << "Arguments: 1. Input file basename" << endl;
    cerr << "           2. Output file basename" << endl;
    cerr << "           3. Nhi  4. rstart  5. rstop" << endl;
    cerr << "           6. m1  7. m2  7. l  9. m" << endl;
    cerr << "           10. Size of l buffer" << endl;
    cerr << "           11. Use horizon flux [Yes(1)/No(0)]" << endl;
    cerr << " Reads a total of Nhi files, ordered from vmax at" << endl;
    cerr << " index 1, to vmin at index Nhi." << endl;
    exit(1);
  }
  char inbase[100], outbase[100], trajname[120], clmname[120];
  FILE *clmfile, *trajfile;
  //
  sprintf(inbase, "%s", argv[1]);
  sprintf(outbase, "%s", argv[2]);
  const int Nhi = atoi(argv[3]);
  const Real rstart = (Real)atof(argv[4]);
  Real rstop = (Real)atof(argv[5]);
  const Real m1 = (Real)atof(argv[6]);
  const Real m2 = (Real)atof(argv[7]);
  const Real M = m1 + m2;
  const Real mu = m1*m2/M;
  const Real eta = mu/M;
  const int l = atoi(argv[8]);
  const int m = atoi(argv[9]);
  const int BUFFER = atoi(argv[10]);
  const int USE_HORIZ = atoi(argv[11]);
  if (USE_HORIZ != 1 && USE_HORIZ != 0) {
    cerr << "11th argument must be yes(1) or no(0)" << endl;
    exit(2);
  }
  //
  CEID ceid(inbase, Nhi, l + BUFFER);
  //
  // Start spiralling on in...
  //
  // bookkeepers
  int DONE_YET = 0;
  int NEAR_END = 0;
  //
  // Initial conditions
  //
  Real t = 0., phi = 0.;
  Real rad = rstart;
  Complex ZIlm, CIlm;
  Real EdotI_lm, EdotH_lm, EdotTot_lm;
  Real hlmplus, hlmcross;
  Real dt = 666.;
  int lines = 0, outfile = 1;
  sprintf(clmname, "%s_M_%1.0e_mu_%3.0e_r0_%2.2lf.h%d%d_%d",
	  outbase, M, mu, rstart, l, m, outfile);
  if (!(clmfile = fopen(clmname, "w"))) {
    cerr << "Error opening " << clmname << endl;
    exit(3);
  }
  fprintf(clmfile, "#t/M  radius/M  orbitalphi  Re(hlm)  Im(hlm)\n");
  ceid.Spline_CIlm(l, m);
  while (!DONE_YET) {
    ceid.Get_rdot_omega(rad);
    ceid.Get_fluxes(rad);
    ceid.Get_CIlm(rad, phi, l, m, CIlm);
    hlmplus = -2.*CIlm.real()/(m*m*ceid.Omega*ceid.Omega);
    hlmcross = 2.*CIlm.imag()/(m*m*ceid.Omega*ceid.Omega);
    //
    // dt is chosen so that we get good resolution on the inspiral.
    //
    if (dt > 1.e-7 && !NEAR_END)
      dt = eta/(10.*fabs(ceid.Omega));
    //
    //  1. (BL time)/M
    //  2. radius/M
    //  3. orbitaphi
    //  4. hlmplus
    //  5. hlmcross
    //
    fprintf(clmfile, "%15.14e %15.14e %15.14e %15.14e %15.14e\n",
    	    t, rad, phi, hlmplus, hlmcross);
    //
    Real rdot, rdot2;
    //
    // 4th order RK step
    //
    if (USE_HORIZ) rdot = ceid.rdot;
    else rdot = ceid.rdot_noH;
    int DONE_RK_STEP = 0;
    Real kr1, kr2, kr3, kr4;
    Real kphi1, kphi2, kphi3, kphi4;
    if (rstop < ceid.rmin) rstop = ceid.rmin;
    while (!DONE_RK_STEP) {
      kr1 = rdot*dt;
      kphi1 = ceid.Omega*dt/eta;
      ceid.Get_rdot_omega(rad + 0.5*kr1);
      if (USE_HORIZ) rdot2 = ceid.rdot;
      else rdot2 = ceid.rdot_noH;
      kr2 = rdot2*dt;
      kphi2 = ceid.Omega*dt/eta;
      ceid.Get_rdot_omega(rad + 0.5*kr2);
      if (USE_HORIZ) rdot2 = ceid.rdot;
      else rdot2 = ceid.rdot_noH;
      kr3 = rdot2*dt;
      kphi3 = ceid.Omega*dt/eta;
      ceid.Get_rdot_omega(rad + kr3);
      if (USE_HORIZ) rdot2 = ceid.rdot;
      else rdot2 = ceid.rdot_noH;
      kr4 = rdot2*dt;
      kphi4 = ceid.Omega*dt/eta;
      //
      if (rad + kr1/6. + kr2/3. + kr3/3. + kr4/6. < rstop ||
	  kr1/6. + kr2/3. + kr3/3. + kr4/6. > 0. ||
	  isnan(rad + kr1/6. + kr2/3. + kr3/3. + kr4/6.)) {
	// Step goes inside rmin: Make dt smaller
	dt /= 2.;
	NEAR_END = 1;
      } else {
	// Step is OK
	DONE_RK_STEP = 1;
      }
    }
    if (dt < 1.e-7) DONE_YET = 1;
    t += dt/eta;
    phi += kphi1/6. + kphi2/3. + kphi3/3. + kphi4/6.;
    rad += kr1/6. + kr2/3. + kr3/3. + kr4/6.;
    lines++;
    if (!(lines%1000000)) {
      fclose(clmfile);
      outfile++;
      sprintf(clmname, "%s_M_%1.0e_mu_%3.0e_r0_%2.2lf.h%d%d_%d",
	      outbase, M, mu, rstart, l, m, outfile);
      if (!(clmfile = fopen(clmname, "w"))) {
	cerr << "Error opening " << clmname << endl;
	exit(3);
      }
    }
  }
  fclose(clmfile);
}

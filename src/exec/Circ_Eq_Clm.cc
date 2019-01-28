//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Clm.cc,v 1.17 2018/08/03 00:21:18 sahughes Exp $
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
  if(argc != 12) {
    cerr << "Arguments: 1. Input file basename" << endl;
    cerr << "           2. Output file basename" << endl;
    cerr << "           3. rmin  4. rmax  5. dr  6. rstart" << endl;
    cerr << "           7. M  8. mu  9. l  10. m" << endl;
    cerr << "           11. Use horizon flux [Yes(1)/No(0)]" << endl;
    exit(1);
  }
  char inbase[40], outbase[40], trajname[40], clmname[80];
  FILE *clmfile, *trajfile;
  //
  sprintf(inbase, "%s", argv[1]);
  sprintf(outbase, "%s", argv[2]);
  const Real rmin = (Real)atof(argv[3]);
  const Real rmax = (Real)atof(argv[4]);
  const Real dr = (Real)atof(argv[5]);
  const Real rstart = (Real)atof(argv[6]);
  const Real M = (Real)atof(argv[7]);
  const Real mu = (Real)atof(argv[8]);
  const Real eta = mu/M;
  const int l = atoi(argv[9]);
  const int m = atoi(argv[10]);
  const int USE_HORIZ = atoi(argv[11]);
  if (USE_HORIZ != 1 && USE_HORIZ != 0) {
    cerr << "11th argument must be yes(1) or no(0)" << endl;
    exit(2);
  }
  //
  CEID ceid(inbase, rmin, rmax, dr, l);
  //
  // Start spiralling on in...
  //
  // bookkeeper
  int DONE_YET = 0;
  //
  // Initial conditions
  //
  Real t = 0., phi = 0.;
  Real rad = rstart;
  Complex Clm;
  Real hlmplus, hlmcross;
  Real dt = 666.;
  int lines = 0, outfile = 1;
  sprintf(clmname, "%s_M_%1.0e_mu_%2.0lf_r0_%2.2lf.h%d%d_%d",
	  outbase, M, mu, rstart, l, m, outfile);
  if (!(clmfile = fopen(clmname, "w"))) {
    cerr << "Error opening " << clmname << endl;
    exit(3);
  }
  sprintf(trajname, "%s_M_%1.0e_mu_%2.0lf_r0_%2.2lf.traj_%d",
	  outbase, M, mu, rstart, outfile);
  if (!(trajfile = fopen(trajname, "w"))) {
    cerr << "Error opening " << trajname << endl;
    exit(4);
  }
  ceid.Spline_CIlm(l, m);
  while (!DONE_YET) {
    ceid.Get_rdot_omega(rad);
    ceid.Get_CIlm(rad, phi, l, m, Clm);
    hlmplus = -2.*Clm.real()/(m*m*ceid.Omega*ceid.Omega);
    hlmcross = 2.*Clm.imag()/(m*m*ceid.Omega*ceid.Omega);
    //
    // dt is chosen so that we get good resolution on the inspiral.
    //
    if (dt > 1.e-7)
      dt = eta/(10.*fabs(ceid.Omega));
    //
    //  1. (BL time)/M
    //  2. hlmplus
    //  3. hlmcross
    //  4. Clm_real()
    //  5. Clm_imag()
    //
    fprintf(clmfile, "%15.14e %15.14e %15.14e %15.14e %15.14e\n",
	    t, hlmplus, hlmcross, eta*Clm.real(), eta*Clm.imag());
    //
    //  1. (BL time)/M
    //  2. rad/M
    //  3. phi
    //  4. M omega
    //
    fprintf(trajfile, "%15.14e %15.14e %15.14e %15.14e\n",
	    t, rad, phi, ceid.Omega);
    //
    // Ending condition: If dt gets way tiny, kill the loop.  We
    // decrease dt if the next step will go inside the minimum
    // radius.
    //
    if (dt > 1.e-7) {
      Real dr_tmp, rdot;
      if (USE_HORIZ) rdot = ceid.rdot;
      else rdot = ceid.rdot_noH;
      dr_tmp = rdot*dt;
      while (rad + dr_tmp < ceid.rmin && dt > 1.e-8) {
	dt /= 2.;
	dr_tmp = rdot*dt;
      }
      //
      // Update the coordinates.
      //
      t += dt/eta;
      phi += ceid.Omega*dt/eta;
      rad += dr_tmp;
      lines++;
      if (!(lines%1000000)) {
	fclose(clmfile);
	fclose(trajfile);
	outfile++;
	sprintf(clmname, "%s_M_%1.0e_mu_%2.0lf_r0_%2.2lf.h%d%d_%d",
		outbase, M, mu, rstart, l, m, outfile);
	if (!(clmfile = fopen(clmname, "w"))) {
	  cerr << "Error opening " << clmname << endl;
	  exit(3);
	}
	sprintf(trajname, "%s_M_%1.0e_mu_%2.0lf_r0_%2.2lf.traj_%d",
		outbase, M, mu, rstart, outfile);
	if (!(trajfile = fopen(trajname, "w"))) {
	  cerr << "Error opening " << trajname << endl;
	  exit(4);
	}
      }
    } else {
      DONE_YET = 1;
    }
  }
  fclose(clmfile);
}

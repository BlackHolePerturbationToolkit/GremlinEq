//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Wave.cc,v 1.15 2018/08/02 19:32:01 sahughes Exp $
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
  if(argc != 10) {
    cerr << "Arguments: 1. data file basename" << endl;
    cerr << "           2. rmin  3. rmax  4. dr  5. rstart" << endl;
    cerr << "           6. lmax  7. M  8. mu  9. cos(th_view)" << endl;
    exit(1);
  }
  char basename[40], wavename[40];
  FILE *wavefile;
  //
  sprintf(basename, "%s", argv[1]);
  const Real rmin = (Real)atof(argv[2]);
  const Real rmax = (Real)atof(argv[3]);
  const Real dr = (Real)atof(argv[4]);
  const Real rstart = (Real)atof(argv[5]);
  const int lmax = atoi(argv[6]);
  const Real M = (Real)atof(argv[7]);
  const Real mu = (Real)atof(argv[8]);
  const Real eta = mu/M;
  Real costheta_view = (Real)atof(argv[9]), phi_view = 0.;
  if (costheta_view > 0.9999)
    costheta_view = 0.9999;
  if (costheta_view < -0.9999)
    costheta_view = -0.9999;
  //
  CEID ceid(basename, rmin, rmax, dr, lmax);
  ceid.Spheroids(costheta_view);
  //
  // Start spiralling on in...
  //
  sprintf(wavename, "%s.wave", basename);
  if (!(wavefile = fopen(wavename, "w"))) {
    cerr << "Error opening " << wavename << endl;
    exit(3);
  }
  // bookkeeper
  int DONE_YET = 0;
  //
  // dt is chosen so that we get good resolution on the inspiral.
  //
  Real dt = eta/(10.*ceid.Omega_max);
  //
  // Initial conditions
  //
  Real t = 0., phi = 0.;
  Real rad = rstart;
  while (!DONE_YET) {
    ceid.Get_rdot_omega(rad);
    ceid.Get_wave(rad, phi, phi_view);
    //
    //  1. (BL time)/M
    //  2. (radius)/M
    //  3. phi
    //  4. (radius)cos(phi)/M
    //  5. (radius)sin(phi)/M
    //  6. hplus
    //  7. hcross
    //
    fprintf(wavefile, "%e %lf %lf %lf %lf %lf %lf\n",
	    t, rad, phi, rad*cos(phi), rad*sin(phi),
	    ceid.hp, ceid.hc);
    //
    // Ending condition: If dt gets way tiny, kill the loop.  We
    // decrease dt if the next step will go inside the minimum
    // radius.
    //
    if (dt > 1.e-7) {
      Real dr_tmp = ceid.rdot*dt;
      while (rad + dr_tmp < ceid.rmin && dt > 1.e-8) {
	dt /= 2.;
	dr_tmp = ceid.rdot*dt;
      }
      //
      // Update the coordinates.
      //
      t += dt/eta;
      phi += ceid.Omega*dt/eta;
      rad += dr_tmp;
    } else {
      DONE_YET = 1;
    }
  }
  fclose(wavefile);
}

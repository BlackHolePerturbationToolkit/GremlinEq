//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Seq2.cc,v 1.17 2018/08/02 16:27:33 sahughes Exp $
//
//---------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Globals.h"
#include "CETD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  char basename[50];
  if(argc != 10) {
    cerr << "Arguments: 1. v_min  2. a" << endl;
    cerr << "           3. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           4. epsilon  5. lmax_min" << endl;
    cerr << "           6. N_tot  7. Nstart  8. Nstop" << endl;
    cerr << "           9. output file basename" << endl;
    cerr << endl;
    cerr << "v_max will be taken to the lesser of v for the orbit" << endl;
    cerr << "at the ISCO, or V_CAP." << endl;
    cerr << endl;
    cerr << "At each orbit, we compute until the sum over l has" << endl;
    cerr << "converged to an accuracy of epsilon, or until l" << endl;
    cerr << "reaches lmax_min, whichever comes last."<< endl;
    cerr << endl;
    cerr << "A total of N_tot orbits will be computed, evenly" << endl;
    cerr << "stepping in v from v_min to v_max.  The count" << endl;
    cerr << "starts at Nstart and ends at Nstop, so that runs" << endl;
    cerr << "can be setup to run in parallel." << endl;
    exit(1);
  }
  const Real v_min = (Real)atof(argv[1]);
  const Real a = (Real)atof(argv[2]);
  const int proret = atoi(argv[3]);
  Real r_isco, v_max;
  if (proret) {
    r_isco = Kerr::isco_pro(a);
    v_max = pow(Kerr::Omega_phi_eqpro(r_isco, a), 1./3.);
  } else {
    r_isco = Kerr::isco_ret(a);
    v_max = pow(fabs(Kerr::Omega_phi_eqret(r_isco, a)), 1./3.);
  }
  const Real EPSILON  = (Real)atof(argv[4]);
  const int lmax_min = atoi(argv[5]);
  const int N_tot = atoi(argv[6]);
  const int N_start = atoi(argv[7]);
  const int N_stop = atoi(argv[8]);
  if (N_start > N_tot) {
    cerr << "N_start > N_tot ... that ain't right." << endl;
    exit(2);
  }
  int i;
  Real v, r;
  const Real dv = (v_max - v_min)/((Real)(N_tot - 1));
  for (i = N_start; i <= N_stop; i++) {
    sprintf(basename, "%s_%d", argv[9], i);
    //
    // i = 1 corresponds to isco: Need to map low i to low r.
    //
    v = v_max - ((Real)(i-1))*dv;
    //
    if (proret)
      r = pow(1./(v*v*v) - a, 2./3.);
    else
      r = pow(1./(v*v*v) + a, 2./3.);
    //
    CETD cetd(proret, r, a, basename);
    //
    cetd.Driver(EPSILON, 70, lmax_min);
  }
}

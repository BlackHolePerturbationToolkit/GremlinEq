//---------------------------------------------------------------------------
//
// $Id: CETD.cc,v 1.29 2019/01/26 14:38:43 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial Teukolsky driver
// Scott Hughes, 26 July 2003
//
#include <cmath>
#include "Globals.h"
#include "CETD.h"
#include "CEDR.h"
#include "GKG.h"
#include "CEKG.h"
#include "SWSH.h"
#include "FT.h"
#include "CEKR.h"
#include "RRGW.h"

CETD::CETD(const int orbitsense, const Real rad, const Real spin, char outbase[]) :
  proret(orbitsense), r(rad), a(spin)
{
  rrgw = new RRGW;
  cekg = new CEKG(proret, r, a);
  //
  // Initialize the range of harmonics to something ridiculous so that in
  // DoHarmonic they get set correctly.
  //
  lmin = 10000; lmax = -1; mmin = 10000; mmax = -10000;
  //
  // Set up the output file.  "params" refers to the geodesic parameters
  // that we save: r, a, E, Lz, Om_phi.
  //
  hsize_t dim[1]; dim[0] = 5;
  Real params[5] = {r, a, cekg->E, cekg->Lz, cekg->Om_phi};
 
  //
  // Check if a file already exists with the same spacetime/orbit parameters. 
  // If so write to this file. Otherwise create a new file
  //
  sprintf(outname, "%s.h5", outbase);
  CEDR cedr(outname);
  if(cedr.exists) {
    hdffile = cedr.infile;
    
    Real file_params[5];
    H5LTread_dataset_double(cedr.infile, "/params/", file_params);
    for(int n = 0; n < 5; n++){
      if(file_params[n] != params[n]) {
	cout << "Output file " << outname << " exists but has different spacetime/orbit parameters. Please choose a different filename." << endl;
	exit(0);
      }
    }
  } else {
    hdffile = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group_id = H5Gcreate2(hdffile, "/modes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5LTmake_dataset_double(hdffile, "/params", 1, dim, params);
  }
  
}

CETD::~CETD()
{
  //
  // Before deleting various objects and closing the file, write out the
  // index range.
  //
  //hsize_t dim[1]; dim[0] = 4;

  //H5LTmake_dataset_int(hdffile, "/indexrange", 1, dim, indexrange);
  //
  delete rrgw;
  delete cekg;
   H5Fclose(hdffile);
}

void CETD::Driver(const int lmax)
{
  for (l = 2; l <= lmax; l++) {
    //
    // Only do positive m; can get negative m by symmetry.
    // Note: Go backwards from m = l, since those are the
    // largest m multipoles.
    //
    for(m = l; m >= 1; m--) {
      DoHarmonic(l, m);
    }
  }
}

void CETD::Driver(const Real EPS_L, const int lmax_min)
{
  int dE_l_small = 0;
  int DONE = 0;
  Real dE_l, dE_l_max = 0;

  l = 2;
  do { // l-loop
    dE_l = 0.;
    //
    // Only do positive m; can get negative m by symmetry.
    // Note: Go backwards from m = l, since those are the
    // largest m multipoles.
    //
    for(m = l; m >= 1; m--) {
      DoHarmonic(l, m);
      dE_l += 2.*(EdotH + EdotI);
    }
    if (dE_l > dE_l_max) dE_l_max = dE_l;
    //
    // The l-loop stops when it's had 3 iterations in a row giving a
    // total energy flux EPS_L smaller than the maximum energy flux
    //
    if (dE_l < EPS_L * dE_l_max)
      dE_l_small++;
    else
      dE_l_small = 0;
    l++;
    if (dE_l_small > 2 && l > lmax_min)
      DONE = 1;
  } while(!DONE);
}

void CETD::Driver(const Real EPS_L, const int lmax, const int lmax_min)
{
  int dE_l_small = 0;
  int DONE = 0;
  Real dE_l, dE_l_max = 0;

  l = 2; 
  do { // l-loop
    dE_l = 0.;
    //
    // Only do positive m; can get negative m by symmetry.
    // Note: Go backwards from m = l, since those are the
    // largest m multipoles.
    //
    for(m = l; m >= 1; m--) {
      DoHarmonic(l, m);
      dE_l += 2.*(EdotH + EdotI);
    }
    if (dE_l > dE_l_max) dE_l_max = dE_l;
    //
    // The l-loop stops when it's had 3 iterations in a row giving a
    // total energy flux EPS_L smaller than the maximum energy flux
    //
    if (dE_l < EPS_L * dE_l_max)
      dE_l_small++;
    else
      dE_l_small = 0;
    l++;
    if (dE_l_small > 2 && l > lmax_min)
      DONE = 1;
    if (l > lmax)
      DONE = 1;
  } while(!DONE);
}

void CETD::DoHarmonic(const int l, const int m)
{
  Real wm = m*cekg->Om_phi;
  Real pm = wm - m*a/(2.*Kerr::rplus(a));
  SWSH swsh(-2, l, m, a*wm);
  //
  // If swsh fails to compute an eigenvalue, it swsh returns -123456
  // for the eigenvalue E.  In this case, just skip the next computations:
  // This harmonic is garbage.  Otherwise, keep computing.
  // (NOTE: This is not true anymore.  Need to look into gsl error handling
  // and recover this functionality.)
  //
  Complex Rin, dRin, Rup, dRup;
  Real LzdotI, LzdotH;
  //
  // Note: Qdot is trivial for the equatorial case, but we compute it
  // (and find zero) in order to recycle code that was originally
  // written for the circular inclined case.
  //
  Real QdotI, QdotH, rdotI, rdotH;
  if ((fabs(swsh.E + 123456.) < 1.e-14)) {
    ZI = 0.; ZH = 0.;
    Rin = 0.; dRin = 0.; Rup = 0.; dRup = 0.;
  } else {
    FT ft(l, m, r, a, wm, swsh.lambda, 3.e-14);
    CEKR cekr(&swsh, &ft, cekg);
    //
    if (ft.Accuracy_in() < 5.e-5) {
      ZI = cekr.ZI;
      Rin = ft.TeukRin()/ft.b_trans;
      dRin = ft.dr_TeukRin()/ft.b_trans;
    } else {
      ZI = 0.;
      Rin = 0.;
      dRin = 0.; 
    }
    if (ft.Accuracy_up() < 5.e-5) {
      ZH = cekr.ZH;
      Rup = ft.TeukRup()/ft.c_trans;
      dRup = ft.dr_TeukRup()/ft.c_trans;
    } else {
      ZH = 0.;
      Rup = 0.;
      dRup = 0.;
    }
    //
    // Further sanity check.  Not sure what thresholds to choose,
    // so picking numbers kind of at random.
    //
    if (m + 5 < l && abs(ZH) > 5.)
      ZH = 0.;
    if (m + 5 < l && abs(ZI) > 1.)
      ZI = 0.;
    if (l < lmin) lmin = l;
    if (l > lmax) lmax = l;
    if (m < mmin) mmin = m;
    if (m > mmax) mmax = m;
  }
  rrgw->Flux_Horizon(a, m, swsh.lambda, wm, pm, ZH, EdotH, LzdotH);
  rrgw->Qdotrdot(r, a, 0.0, cekg->E, cekg->Lz, EdotH, LzdotH, QdotH, rdotH);
  rrgw->Flux_Infinity(a, m, swsh.lambda, wm, pm, ZI, EdotI, LzdotI);
  rrgw->Qdotrdot(r, a, 0.0, cekg->E, cekg->Lz, EdotI, LzdotI, QdotI, rdotI);
  //
  // Set up the dataset
  //
  char dataset[50];
  sprintf(dataset, "modes/l%dm%d", l, m);
  //
  // Format the data to go into the dataset
  //
  hsize_t dim[1]; dim[0] = 19;
  Real modedata[19] = {swsh.lambda,
		       ZI.real(), ZI.imag(), ZH.real(), ZH.imag(),
		       Rin.real(), Rin.imag(), dRin.real(), dRin.imag(),
		       Rup.real(), Rup.imag(), dRup.real(), dRup.imag(),
		       EdotI, EdotH, LzdotI, LzdotH, rdotI, rdotH};
  
  
	if(H5Lexists(hdffile, dataset, H5P_DEFAULT)){
		hid_t dataset_id = H5Dopen2(hdffile, dataset, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, modedata);
		H5Dclose(dataset_id);
	}else{
        H5LTmake_dataset_double(hdffile, dataset, 1, dim, modedata);
    }
	
	if(H5Lexists(hdffile, "/indexrange", H5P_DEFAULT)){
		int indexrange[4] = {2, 10, 1, 10};
		hid_t dataset_id = H5Dopen2(hdffile, "/indexrange", H5P_DEFAULT);
		
		H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, indexrange);
		if(l < indexrange[0]) indexrange[0] = l;
		if(l > indexrange[1]) indexrange[1] = l;
		if(m < indexrange[2]) indexrange[2] = m;
		if(m > indexrange[3]) indexrange[3] = m;
			
		H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, indexrange);
		H5Dclose(dataset_id);
		
	}else{
      int indexrange[4] = {l,l,m,m}; 		
  	  dim[0] = 4;
  	  H5LTmake_dataset_int(hdffile, "/indexrange", 1, dim, indexrange);
		
	}
	
  	//FIXME: move indexrange stuff here
  
  
  //
  // Below output is just to follow along on long runs.  Comment out if you don't
  // care to read it.
  //
  Real v = pow(fabs(cekg->Om_phi), 1./3.);
  Real EdotNorm = (32./5.)*pow(v, 10.);
  if (ZI == 0. || ZH == 0.) {
    cout << r << " " << v << " " << l << " " << m << " " << wm << " "
   	 << EdotH/EdotNorm << " " << EdotI/EdotNorm << " "
   	 << "MODE ZEROED" << endl;
  } else {
    cout << r << " " << v << " " << l << " " << m << " " << wm << " "
   	 << EdotH/EdotNorm << " " << EdotI/EdotNorm << endl;
  }
}

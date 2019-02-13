//---------------------------------------------------------------------------
//
// $Id: CEDR.cc,v 1.2 2018/08/04 16:00:53 sahughes Exp $
//
//---------------------------------------------------------------------------
//
//
// HDF5 format data reader for circular equatorial case
//
#include <cmath>
#include <sys/stat.h>
#include <hdf5_hl.h>
#include "Globals.h"
#include "CEDR.h"

CEDR::CEDR(char inname[])
{
  if(fileExists(inname)){
	exists = true;
    infile = H5Fopen(inname, H5F_ACC_RDWR, H5P_DEFAULT);
    //
    // Orbit parameters
    //
    Real params[5];
    H5LTread_dataset_double(infile, "/params/", params);
    r = params[0];
    a = params[1];
    E = params[2];
    Lz = params[3];
    Om_phi = params[4];
    //
    // Range over which the indices span
    //
    int indexrange[4];
    H5LTread_dataset_int(infile, "/indexrange", indexrange);
    lmin = indexrange[0];
    lmax = indexrange[1];
    mmin = indexrange[2];
    mmax = indexrange[3];
	
  }else{
	exists = false;
  }
}

CEDR::~CEDR()
{
	//H5Fclose(infile);
}

// Check if a file exists
bool CEDR::fileExists(const std::string& filename){
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1) return true;
    return false;
}

int CEDR::ReadData(const int l, const int m)
{
  char modeset[50];
  Real modedata[19];
  if (l >= lmin && l <= lmax &&
      m >= mmin && m <= mmax && abs(m) <= l) {
    sprintf(modeset, "/modes/l%dm%d", l, m);
    if (H5Lexists(infile, modeset, H5P_DEFAULT)) {
      H5LTread_dataset_double(infile, modeset, modedata);
      //
      lambda = modedata[0];
      ZI = Complex(modedata[1], modedata[2]);
      ZH = Complex(modedata[3], modedata[4]);
      Rin = Complex(modedata[5], modedata[6]);
      dRin = Complex(modedata[7], modedata[8]);
      Rup = Complex(modedata[9], modedata[10]);
      dRup = Complex(modedata[11], modedata[12]);
      EdotI = modedata[13];
      EdotH = modedata[14];
      LzdotI = modedata[15];
      LzdotH = modedata[16];
      rdotI = modedata[17];
      rdotH = modedata[18];
      return(1);
    } else return(0);
  } else return(0);
}

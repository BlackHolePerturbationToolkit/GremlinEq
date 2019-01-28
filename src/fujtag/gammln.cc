//---------------------------------------------------------------------------
//
// $Id: gammln.cc,v 1.1 2013/05/17 15:39:53 sahughes Exp $
//
//---------------------------------------------------------------------------

/* gammln.cc
 * the logarithm of the gamma functions, for complex values
 * Based on Numerical Recipies
 */

#ifndef FT_STANDALONE
#  include "FT.h"
#else
#  include "funcs.h"
#  define I std::complex<double>(0,1)
#endif

std::complex<double> FT_PRFX gammln(const std::complex<double> x) {
  int j;
  std::complex<double> tmp, y, ser;
  static const double coef[14] = 
    { 57.1562356658629235,    -59.5979603554754912,     14.1360979747417471,
     -0.491913816097620199,    3.39946499848118887e-5,  4.65236289270485756e-5,
     -9.83744753048795646e-5,  1.58088703224912494e-4, -2.10264441724104883e-4,
      2.17439618115212643e-4, -1.64318106536763890e-4,  8.44182239838527433e-5,
     -2.61908384015814087e-5,  3.68991826595316234e-6};
  if(real(x) <= 0) return log(M_PI) - sinln(M_PI * x) - gammln(1.0 - x);
  y = x;
  tmp = x + 5.24218750000000000;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for(j=0;j<14;j++) ser += coef[j]/(y += 1.0);
  return tmp+log(2.5066282746310005*ser/x);
}

std::complex<double> FT_PRFX sinln(const std::complex<double> x) {
  return - I * M_PI_2 - M_LN2 +
    ((imag(x) > 0) ?
     - I * x + log(exp(2.*I*x)-1.) :
     + I * x + log(1.-exp(-2.*I*x)));
}

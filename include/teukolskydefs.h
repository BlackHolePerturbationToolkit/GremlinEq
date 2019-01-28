//---------------------------------------------------------------------------
//
// $Id: teukolskydefs.h,v 1.2 2014/04/08 19:22:54 sahughes Exp $
//
//---------------------------------------------------------------------------


/* teukolskydefs.h
 * general definitions for the Teukolsky project
 */

#ifndef TEUKOLSKYDEFS_H_SEEN
#define TEUKOLSKYDEFS_H_SEEN

#include <cfloat>
#include <complex>
#include "complex_ops.h"

extern "C" int isnan(double);

/* The precision to use
 * 0 float
 * 1 double
 * 2 long double
 */
#ifndef REAL_TYPE
#  define REAL_TYPE 1
#endif

#if REAL_TYPE == 0
typedef float Real;
#define REAL_MANT_DIG FLT_MANT_DIG
#if DBL_MANT_DIG > FLT_MANT_DIG
typedef double LongReal;
#define LREAL_MANT_DIG DBL_MANT_DIG
#else
typedef long double LongReal;
#define LREAL_MANT_DIG LDBL_MANT_DIG
#endif

#define fsRe "hf"
#define REAL_EPSILON FLT_EPSILON
#define PREC_FUNCSUFFIX f

#elif REAL_TYPE == 1
typedef double Real;
#define REAL_MANT_DIG DBL_MANT_DIG
typedef long double LongReal;
#define LREAL_MANT_DIG LDBL_MANT_DIG

#define fsRe "f"
#define REAL_EPSILON DBL_EPSILON
#define PREC_FUNCSUFFIX

#elif REAL_TYPE == 2
typedef long double Real;
#define REAL_MANT_DIG LDBL_MANT_DIG
typedef long double LongReal;
#define LREAL_MANT_DIG LDBL_MANT_DIG

#define fsRe "Lf"
#define REAL_EPSILON LDBL_EPSILON
#define PREC_FUNCSUFFIX l

#else
#error "Bad REAL_TYPE"
#endif /* REAL_TYPE */

typedef std::complex<Real> Complex;
typedef std::complex<LongReal> LongComplex;
#define I (Complex(0,1))

/* precision dependent functions */
#define MAKE_PRECFUNC(func) FUNC_CONCAT(func,PREC_FUNCSUFFIX)
#define FUNC_CONCAT(a,b) FUNC_CONCAT_INNER(a,b)
#define FUNC_CONCAT_INNER(a,b) a ## b

#define Fmax MAKE_PRECFUNC(fmax)
#define Fmin MAKE_PRECFUNC(fmin)
#define Fabs MAKE_PRECFUNC(fabs)
#define Sqrt MAKE_PRECFUNC(sqrt)

/* useful utilities */
#define CNAN Complex(NAN,NAN)

#define Cisnan(x) (isnan(real(x)) || isnan(imag(x)))
#define Cmaxpart(x) (Fmax(Fabs(real(x)),Fabs(imag(x))))

#endif /* TEUKOLSKYDEFS_H_SEEN */

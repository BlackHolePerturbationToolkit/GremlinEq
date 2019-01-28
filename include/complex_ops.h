//---------------------------------------------------------------------------
//
// $Id: complex_ops.h,v 1.2 2014/08/01 17:17:18 sahughes Exp $
//
//---------------------------------------------------------------------------

/* complex_ops.h
 * add some more operations on complex numbers
 */

#ifndef COMPLEX_OPS_H_SEEN
#define COMPLEX_OPS_H_SEEN

#include <complex>
namespace std {

  /* math ops with real types (or anything convertable to a float) */

  /* general cases */

  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator+(const complex<_Tp>& __x, const _S& __y) {
    complex<_Tp> __r = __x;
    __r += __y;
    return __r;
  }

  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator+(const _S& __x, const complex<_Tp>& __y) {
    complex<_Tp> __r = __y;
    __r += __x;
    return __r;
  }


  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator*(const complex<_Tp>& __x, const _S& __y) {
    complex<_Tp> __r = __x;
    __r *= __y;
    return __r;
  }

  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator*(const _S& __x, const complex<_Tp>& __y) {
    complex<_Tp> __r = __y;
    __r *= __x;
    return __r;
  }


  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator-(const complex<_Tp>& __x, const _S& __y) {
    complex<_Tp> __r = __x;
    __r -= __y;
    return __r;
  }

  /* Edited by SAH 1 Aug 2014, new compilers barf on old code */
  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator-(const _S& __x, const complex<_Tp>& __y) {
    complex<_Tp> __r(__x - __y.real(), - __y.imag());
    /* complex<_Tp> __r(__x, - __y.imag()); */
    /* __r.real() -= __y.real(); */
    return __r;
  }


  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator/(const complex<_Tp>& __x, const _S& __y) {
    complex<_Tp> __r = __x;
    __r /= __y;
    return __r;
  }

  template<typename _Tp, typename _S>
    inline complex<_Tp>
    operator/(const _S& __x, const complex<_Tp>& __y) {
    complex<_Tp> __r = __x;
    __r /= __y;
    return __r;
  }

  /* special cases for operations between two complexes of different widths
   * (only double + long double implemented) */

  inline complex<long double>
    operator+(const complex<long double>& __x, const complex<double>& __y) {
    complex<long double> __r = __x;
    __r += __y;
    return __r;
  }
  inline complex<long double>
    operator+(const complex<double>& __x, const complex<long double>& __y) {
    complex<long double> __r = __y;
    __r += __x;
    return __r;
  }

  inline complex<long double>
    operator-(const complex<long double>& __x, const complex<double>& __y) {
    complex<long double> __r = __x;
    __r -= __y;
    return __r;
  }
  inline complex<long double>
    operator-(const complex<double>& __x, const complex<long double>& __y) {
    complex<long double> __r = - __y;
    __r += __x;
    return __r;
  }

  inline complex<long double>
    operator*(const complex<long double>& __x, const complex<double>& __y) {
    complex<long double> __r = __x;
    __r *= __y;
    return __r;
  }
  inline complex<long double>
    operator*(const complex<double>& __x, const complex<long double>& __y) {
    complex<long double> __r = __y;
    __r *= __x;
    return __r;
  }

  inline complex<long double>
    operator/(const complex<long double>& __x, const complex<double>& __y) {
    complex<long double> __r = __x;
    __r /= __y;
    return __r;
  }
  inline complex<long double>
    operator/(const complex<double>& __x, const complex<long double>& __y) {
    complex<long double> __r = 1/__y;
    __r *= __x;
    return __r;
  }
}

#endif /* !COMPLEX_OPS_H_SEEN */

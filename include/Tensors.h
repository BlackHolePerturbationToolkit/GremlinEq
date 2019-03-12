//---------------------------------------------------------------------------
//
// $Id: Tensors.h,v 1.3 2019/01/24 20:03:19 sahughes Exp $
//
//---------------------------------------------------------------------------
//
// Class Tensor ... a container for defining multi-index array and tensor
// type objects in a numerical-recipes like fashion, but good for any
// type.  Note all the methods are static functions, so you don't need
// create an instance of the class to call them.
//
// Example: To make a 4-index tensor of type Complex, with indices
// running (al, ah; bl, bh; cl, ch; dl, dh):
//
// Complex ****Z;
// Z = Tensor<Complex>::tensor4(al, ah, bl, bh, cl, ch, dl, dh);
//
// To free the memory associated with this:
// Tensor<Complex>::free_tensor4(Z, al, ah, bl, bh, cl, ch, dl, dh);
//

#ifndef _Tensors_H
#define _Tensors_H

#include "Globals.h"

#define OFFSET 0

template<class TypeHere>

//! Tensors Class
/*! This is a container class for memory allocation routines which define multi-index objects with arbitrary index range. They are inspired by Numerical Recipe’s such routines [e.g., vector() and matrix()], but have been implemented in a new way and use C++ templates to make arrays of arbitrary type. */
class Tensor {
 public:
  static TypeHere *vector(const long al, const long ah)
  {
    TypeHere *t;
    long da = ah - al + 1;
    
    t = (TypeHere *)malloc((size_t)((da + OFFSET)*sizeof(TypeHere)));
    
    return t - al + OFFSET;
  }

  static void free_vector(TypeHere *t,
			  const long al, const long ah)
  {
    free(t + al - OFFSET);
  }

  static TypeHere **vectorptr(const long al, const long ah)
  {
    TypeHere **t;
    long da = ah - al + 1;
    
    t = (TypeHere **)malloc((size_t)((da + OFFSET)*sizeof(TypeHere)));
    
    return t - al + OFFSET;
  }
  
  static void free_vectorptr(TypeHere **t,
			     const long al, const long ah)
  {
    free(t + al - OFFSET);
  }
  
  static TypeHere **matrix(const long al, const long ah,
			   const long bl, const long bh)
  {
    TypeHere **t;
    long i, da = ah - al + 1;
    
    t = (TypeHere **)malloc((size_t)(da + OFFSET)*sizeof(TypeHere*));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = vector(bl, bh);
    
    return t;
  }
  
  static void free_matrix(TypeHere **t,
			  const long al, const long ah,
			  const long bl, const long bh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_vector(t[i], bl, bh);
    
    free(t + al - OFFSET);
  }
    
  static TypeHere ***matrixptr(const long al, const long ah,
			       const long bl, const long bh)
  {
    TypeHere ***t;
    long i, da = ah - al + 1;
    
    t = (TypeHere ***)malloc((size_t)(da + OFFSET)*sizeof(TypeHere**));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = vectorptr(bl, bh);
    
    return t;
  }
  
  static void free_matrixptr(TypeHere ***t,
			     const long al, const long ah,
			     const long bl, const long bh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_vectorptr(t[i], bl, bh);
    
    free(t + al - OFFSET);
  }
  
  static TypeHere ***tensor3(const long al, const long ah,
			     const long bl, const long bh,
			     const long cl, const long ch)
  {
    TypeHere ***t;
    long i, da = ah - al + 1;
    
    t = (TypeHere ***)malloc((size_t)(da + OFFSET)*sizeof(TypeHere**));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = matrix(bl, bh, cl, ch);
    
    return t;
  }
  
  static void free_tensor3(TypeHere ***t,
			   const long al, const long ah,
			   const long bl, const long bh,
			   const long cl, const long ch)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_matrix(t[i], bl, bh, cl, ch);
    
    free(t + al - OFFSET);
  }

  static TypeHere ****tensor4(const long al, const long ah,
			      const long bl, const long bh,
			      const long cl, const long ch,
			      const long dl, const long dh)
  {
    TypeHere ****t;
    long i, da = ah - al + 1;
    
    t = (TypeHere ****)malloc((size_t)(da + OFFSET)*sizeof(TypeHere***));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = tensor3(bl, bh, cl, ch, dl, dh);
    
    return t;
  }
  
  static void free_tensor4(TypeHere ****t,
			   const long al, const long ah,
			   const long bl, const long bh,
			   const long cl, const long ch,
			   const long dl, const long dh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_tensor3(t[i], bl, bh, cl, ch, dl, dh);
    
    free(t + al - OFFSET);
  }

  static TypeHere *****tensor5(const long al, const long ah,
			       const long bl, const long bh,
			       const long cl, const long ch,
			       const long dl, const long dh,
			       const long el, const long eh)
  {
    TypeHere *****t;
    long i, da = ah - al + 1;
    
    t = (TypeHere *****)malloc((size_t)(da + OFFSET)*sizeof(TypeHere****));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = tensor4(bl, bh, cl, ch, dl, dh, el, eh);
    
    return t;
  }
  
  static void free_tensor5(TypeHere *****t,
			   const long al, const long ah,
			   const long bl, const long bh,
			   const long cl, const long ch,
			   const long dl, const long dh,
			   const long el, const long eh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_tensor4(t[i], bl, bh, cl, ch, dl, dh, el, eh);
    
    free(t + al - OFFSET);
  }

  static TypeHere ******tensor6(const long al, const long ah,
				const long bl, const long bh,
				const long cl, const long ch,
				const long dl, const long dh,
				const long el, const long eh,
				const long fl, const long fh)
  {
    TypeHere ******t;
    long i, da = ah - al + 1;
    
    t = (TypeHere ******)malloc((size_t)(da + OFFSET)*sizeof(TypeHere*****));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = tensor5(bl, bh, cl, ch, dl, dh, el, eh, fl, fh);
    
    return t;
  }
  
  static void free_tensor6(TypeHere ******t,
			   const long al, const long ah,
			   const long bl, const long bh,
			   const long cl, const long ch,
			   const long dl, const long dh,
			   const long el, const long eh,
			   const long fl, const long fh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_tensor5(t[i], bl, bh, cl, ch, dl, dh, el, eh, fl, fh);
    
    free(t + al - OFFSET);
  }

  static TypeHere *******tensor7(const long al, const long ah,
				 const long bl, const long bh,
				 const long cl, const long ch,
				 const long dl, const long dh,
				 const long el, const long eh,
				 const long fl, const long fh,
				 const long gl, const long gh)
  {
    TypeHere *******t;
    long i, da = ah - al + 1;
    
    t = (TypeHere *******)malloc((size_t)(da + OFFSET)*sizeof(TypeHere******));
    t += OFFSET;
    t -= al;
    
    for (i = al; i <= ah; i++)
      t[i] = tensor6(bl, bh, cl, ch, dl, dh, el, eh, fl, fh, gl, gh);
    
    return t;
  }
  
  static void free_tensor7(TypeHere *******t,
			   const long al, const long ah,
			   const long bl, const long bh,
			   const long cl, const long ch,
			   const long dl, const long dh,
			   const long el, const long eh,
			   const long fl, const long fh,
			   const long gl, const long gh)
  {
    long i;
    for (i = ah; i >= al; i--)
      free_tensor6(t[i], bl, bh, cl, ch, dl, dh, el, eh, fl, fh, gl, gh);
    
    free(t + al - OFFSET);
  }
};

#endif

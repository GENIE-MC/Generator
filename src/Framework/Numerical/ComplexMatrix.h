//____________________________________________________________________________
/*!

  \namespace  genie::utils::math

  \brief      Utilities for ComplexMatrix definitions

  \author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
  University of Liverpool & STFC Rutherford Appleton Lab

  \created    July, 2020

  \cpright    Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COMPLEX_MATRIX_H_
#define _COMPLEX_MATRIX_H_

#include <vector>
#include <array>
#include <complex>

namespace genie {
  namespace utils {

    namespace math  {

      template< size_t N >
      class ComplexArray : public std::array< std::complex<double>, N > {

        const ComplexArray & operator += ( const ComplexArray & ) ;
        const ComplexArray & operator *= ( const std::complex<double> & ) ;

        ComplexArray Conj() const ;
      } ;

      template< size_t N >
      class ComplexMatrix : public std::array< ComplexArray<N>, N > {

        const ComplexMatrix & operator += ( const ComplexMatrix & ) ;
        const ComplexMatrix & operator *= ( const std::complex<double> & ) ;

        ComplexMatrix Conj() const ;

      };

    } // math  namespace
  } // utils namespace
} // genie namespace

#endif // _COMPLEX_MATRIX_H_

//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 */
//____________________________________________________________________________


#include "Framework/Numerical/ComplexMatrix.h"

using namespace genie::utils::math ;

template<size_t N>
const ComplexArray & ComplexArray::operator += ( const ComplexArray & v ) {

  for ( unsigned int i = 0 ; i < size(); ++i ) {
    at(i) += v[i] ;
  }

  return *this ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexArray & ComplexArray::operator *= ( const complex<double> & c ) {

  for ( unsigned int i = 0 ; i < size(); ++i ) {
    at(i) *= c ;
  }

  return *this ;

}

//____________________________________________________________________________
template<size_t>
ComplexArray ComplexArray::Conj() const {

  ComplexArray a;
  for ( unsigned int i = 0 ; i < size(); ++i ) {
    a.[i] = std::conj( at(i) ) ;
  }

  return a ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexMatrix & ComplexMatrix::operator += ( const ComplexMatrix & m ) {

  for ( unsigned int i = 0 ; i < size(); ++i ) {
    at(i) += v[i] ;
  }

  return *this ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexMatrix & operator *= ( const complex<double> & c ) {

  for ( unsigned int i = 0 ; i < size(); ++i ) {
    at(i) *= c ;
  }

  return *this ;

}

//____________________________________________________________________________
template<size_t N >
ComplexMatrix ComplexMatrix::Conj() const {

  ComplexMatrix a;
  for ( unsigned int i = 0 ; i < size(); ++i ) {
    a.[i] = at(i).Conj() ;
  }

  return a ;
}

//____________________________________________________________________________

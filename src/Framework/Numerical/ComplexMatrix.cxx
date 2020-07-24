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

using std::complex ;

template<size_t N>
const ComplexArray<N> & ComplexArray<N>::operator += ( const ComplexArray<N> & v ) {

  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    this->at(i) += v[i] ;
  }

  return *this ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexArray<N> & ComplexArray<N>::operator *= ( const complex<double> & c ) {

  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    this -> at(i) *= c ;
  }

  return *this ;

}

//____________________________________________________________________________
template<size_t N>
ComplexArray<N> ComplexArray<N>::Conj() const {

  ComplexArray<N> a;
  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    a[i] = std::conj( this->at(i) ) ;
  }

  return a ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexMatrix<N> & ComplexMatrix<N>::operator += ( const ComplexMatrix<N> & m ) {

  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    this->at(i) += m[i] ;
  }

  return *this ;
}

//____________________________________________________________________________

template<size_t N>
const ComplexMatrix<N> & ComplexMatrix<N>::operator *= ( const complex<double> & c ) {

  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    this->at(i) *= c ;
  }

  return *this ;
}

//____________________________________________________________________________
template<size_t N >
ComplexMatrix<N> ComplexMatrix<N>::Conj() const {

  ComplexMatrix<N> a;
  for ( unsigned int i = 0 ; i < this->size(); ++i ) {
    a[i] = this->at(i).Conj() ;
  }

  return a ;
}

//____________________________________________________________________________

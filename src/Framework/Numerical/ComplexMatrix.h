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

      public: 
        const ComplexArray & operator += ( const ComplexArray & ) ;
        const ComplexArray & operator *= ( const std::complex<double> & ) ;

        ComplexArray Conj() const ;
      } ;

      template< size_t N >
	class ComplexMatrix : public std::array< ComplexArray<N>, N > {

      public: 
        const ComplexMatrix & operator += ( const ComplexMatrix & ) ;
        const ComplexMatrix & operator *= ( const std::complex<double> & ) ;

        ComplexMatrix Conj() const ;

      };

      // all these classes are templated, so they need the implementation locally
      

      template<size_t N>
	const ComplexArray<N> & ComplexArray<N>::operator += ( const ComplexArray<N> & v ) {

	for ( unsigned int i = 0 ; i < this->size(); ++i ) {
	  this->at(i) += v[i] ;
	}

	return *this ;
      }

      //____________________________________________________________________________                    

      template<size_t N>
	const ComplexArray<N> & ComplexArray<N>::operator *= ( const std::complex<double> & c ) {

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
	const ComplexMatrix<N> & ComplexMatrix<N>::operator *= ( const std::complex<double> & c ) {

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





    } // math  namespace
  } // utils namespace
} // genie namespace

#endif // _COMPLEX_MATRIX_H_

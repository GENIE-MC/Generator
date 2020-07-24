//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda <mroda \at liverpool.ac.uk> 
         University of Liverpool 

 For documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Numerical/GTrace.h"

using namespace genie::utils::math;

//____________________________________________________________________________
const GTrace &GTrace::operator+=( const GTrace & M ) {
  for( unsigned int i=0; i<M.size(); i++ ) {
    for( unsigned int j=0; j<M.size(); j++ ) {
      (*this)[i][j] += M[i][j];
	}
  }


  return *this;
}
//____________________________________________________________________________
const GTrace &GTrace::operator*=( std::complex<double> c) {
  for( unsigned int i=0; i<(*this).size(); i++ ) {
    for( unsigned int j=0; j<(*this).size(); j++ ) {
      (*this)[i][j] *= c;
	}
  }

  return *this;
}
//____________________________________________________________________________
GTrace GTrace::Conj( const GTrace & M ) {
  for( unsigned int i=0; i<M.size(); i++ ) {
    for( unsigned int j=0; j<M.size(); j++ ) {
      (*this)[i][j] = std::conj(M[i][j]);
    }
  }
  return (*this);

}

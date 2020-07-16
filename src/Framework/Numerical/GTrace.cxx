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
  for( int i=0; i<M.size(); i++; ) {
    for( int j=0; j<M.size(); j++; ) {
      (*this)[i][j] += M[i][j];
	}
  }
  return *this;
}
//____________________________________________________________________________
const GTrace &GTrace::operator*=( std::complex<double> c) {
  for( int i=0; i<M.size(); i++; ) {
    for( int j=0; j<M.size(); j++; ) {
      (*this)[i][j] *= c;
	}
  }
  return *this;
}

//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>

#include "Physics/Coherent/XSection/ARWavefunction.h"

namespace genie
{
namespace alvarezruso
{
    
ARWavefunction::ARWavefunction(unsigned int sampling_in, bool debug)
  : debug_(debug),
  sampling_(2*sampling_in),
  wavefunction_(sampling_, std::vector<std::complex<double> >(sampling_, std::complex<double> (0.0,0.0)) )
{
  if(debug_) std::cerr << "WF@ Constructor" << std::endl;
}

ARWavefunction::~ARWavefunction() {}

std::string ARWavefunction::print() const
{
  std::ostringstream oss;
  oss << "{";
  for(unsigned int i = 0; i != sampling_; ++i)
  {
    oss << "[";
    for(unsigned int j = 0; j != sampling_; ++j)
    {
      oss << ((*this)[i][j]);
      if( j != (sampling_ - 1) ) oss << ", ";
    }
    oss << "]";
  }
  oss << "}";
  return oss.str();
}

const std::vector<std::complex<double> >& ARWavefunction::operator[] (unsigned int i) const
{
  return wavefunction_[i];
}

const std::complex<double> & ARWavefunction::operator() (unsigned int i, unsigned int j) const
{
  return wavefunction_[i][j];
}

std::complex<double>  ARWavefunction::get(unsigned int i, unsigned int j) const
{
  return wavefunction_[i][j];
}

void ARWavefunction::set(unsigned int i, unsigned int j, const std::complex<double> & value)
{
  wavefunction_[i][j] = value;
}

unsigned int ARWavefunction::sampling() const  {  return sampling_;  }

} //namespace alvarezruso
} //namespace genie


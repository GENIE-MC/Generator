//____________________________________________________________________________
/*
 Copyright (c) 2023-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/RESVectFFAmplitude.h"

using namespace genie;
using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream & stream, const RESVectFFAmplitude & hamp)
  {
     hamp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
RESVectFFAmplitude::RESVectFFAmplitude()
{
  this->Init();
}
//____________________________________________________________________________
RESVectFFAmplitude::RESVectFFAmplitude(const RESVectFFAmplitude & hamp)
{
  fA12 = hamp.AmplA12();
  fA32  = hamp.AmplA32();
  fS12 = hamp.AmplS12();
}
//____________________________________________________________________________
void RESVectFFAmplitude::Print(ostream & stream) const
{
  stream << endl;
  stream << " A1/2 = " << fA12 << endl;
  stream << " A3/2 = " << fA32  << endl;
  stream << " S1/2 = " << fS12 << endl;
}
//____________________________________________________________________________
void RESVectFFAmplitude::Init(void)
{
  fA12 = 0. ;
  fA32 = 0. ;
  fS12 = 0. ;
}
//____________________________________________________________________________

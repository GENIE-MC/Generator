//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/RSHelicityAmpl.h"

using namespace genie;
using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream & stream, const RSHelicityAmpl & hamp)
  {
     hamp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
RSHelicityAmpl::RSHelicityAmpl()
{
  this->Init();
}
//____________________________________________________________________________
RSHelicityAmpl::RSHelicityAmpl(const RSHelicityAmpl & hamp)
{
  fMinus1 = hamp.AmpMinus1();
  fPlus1  = hamp.AmpPlus1();
  fMinus3 = hamp.AmpMinus3();
  fPlus3  = hamp.AmpPlus3();
  f0Minus = hamp.Amp0Minus();
  f0Plus  = hamp.Amp0Plus();
}
//____________________________________________________________________________
void RSHelicityAmpl::Print(ostream & stream) const
{
  stream << endl;
  stream << " f(-1) = " << fMinus1 << endl;
  stream << " f(+1) = " << fPlus1  << endl;
  stream << " f(-3) = " << fMinus3 << endl;
  stream << " f(+3) = " << fPlus3  << endl;
  stream << " f(0-) = " << f0Minus << endl;
  stream << " f(0+) = " << f0Plus  << endl;
}
//____________________________________________________________________________
void RSHelicityAmpl::Init(void)
{
  fMinus1 = 0.0;
  fPlus1  = 0.0;
  fMinus3 = 0.0;
  fPlus3  = 0.0;
  f0Minus = 0.0;
  f0Plus  = 0.0;
}
//____________________________________________________________________________



//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
         
 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/MKHelicityAmpl.h"

using namespace genie;
using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream & stream, const MKHelicityAmpl & hamp)
  {
     hamp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
MKHelicityAmpl::MKHelicityAmpl()
{
  this->Init();
}
//____________________________________________________________________________
MKHelicityAmpl::MKHelicityAmpl(const MKHelicityAmpl & hamp)
{
  fVMinus1   = hamp.AmpVMinus1();
  fVPlus1    = hamp.AmpVPlus1();
  fVMinus3   = hamp.AmpVMinus3();
  fVPlus3    = hamp.AmpVPlus3();
  fV0LMinus  = hamp.AmpV0LMinus();
  fV0LPlus   = hamp.AmpV0LPlus();
  fV0RMinus  = hamp.AmpV0RMinus();
  fV0RPlus   = hamp.AmpV0RPlus();
  
  fAMinus1   = hamp.AmpAMinus1();
  fAPlus1    = hamp.AmpAPlus1();
  fAMinus3   = hamp.AmpAMinus3();
  fAPlus3    = hamp.AmpAPlus3();
  fA0LMinus  = hamp.AmpA0LMinus();
  fA0LPlus   = hamp.AmpA0LPlus();
  fA0RMinus  = hamp.AmpA0RMinus();
  fA0RPlus   = hamp.AmpA0RPlus();
  
}
//____________________________________________________________________________
void MKHelicityAmpl::Print(ostream & stream) const
{
  stream << endl;
  stream << " fV(-1)   = " << fVMinus1  << endl;
  stream << " fV(+1)   = " << fVPlus1   << endl;
  stream << " fV(-3)   = " << fVMinus3  << endl;
  stream << " fV(+3)   = " << fVPlus3   << endl;
  stream << " fV(0L-)  = " << fV0LMinus << endl;
  stream << " fV(0L+)  = " << fV0LPlus  << endl;
  stream << " fV(0R-)  = " << fV0RMinus << endl;
  stream << " fV(0R+)  = " << fV0RPlus  << endl;
  
  stream << " fA(-1)   = " << fAMinus1  << endl;
  stream << " fA(+1)   = " << fAPlus1   << endl;
  stream << " fA(-3)   = " << fAMinus3  << endl;
  stream << " fA(+3)   = " << fAPlus3   << endl;
  stream << " fA(0L-)  = " << fA0LMinus << endl;
  stream << " fA(0L+)  = " << fA0LPlus  << endl;
  stream << " fA(0R-)  = " << fA0RMinus << endl;
  stream << " fA(0R+)  = " << fA0RPlus  << endl;
  
}
//____________________________________________________________________________
void MKHelicityAmpl::Init(void)
{
  fVMinus1   = 0.0;
  fVPlus1    = 0.0;
  fVMinus3   = 0.0;
  fVPlus3    = 0.0;
  fV0LMinus  = 0.0;
  fV0LPlus   = 0.0;
  fV0RMinus  = 0.0;
  fV0RPlus   = 0.0;
  
  fAMinus1   = 0.0;
  fAPlus1    = 0.0;
  fAMinus3   = 0.0;
  fAPlus3    = 0.0;
  fA0LMinus  = 0.0;
  fA0LPlus   = 0.0;
  fA0RMinus  = 0.0;
  fA0RPlus   = 0.0;
  
}
//____________________________________________________________________________



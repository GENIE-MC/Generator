//____________________________________________________________________________
/*
Dark Matter Calculations of |f|^2 terms for the cross section, based on the
procedure of FKR.

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplDM.h"

using namespace genie;
using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream & stream, const RSHelicityAmplDM & hamp)
  {
     hamp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
RSHelicityAmplDM::RSHelicityAmplDM()
{
  this->Init();
}
//____________________________________________________________________________
RSHelicityAmplDM::RSHelicityAmplDM(const RSHelicityAmplDM & hamp)
{
  fMinus1 = hamp.AmpMinus1();
  fPlus1  = hamp.AmpPlus1();
  fMinus3 = hamp.AmpMinus3();
  fPlus3  = hamp.AmpPlus3();
  f0Minus = hamp.Amp0Minus();
  f0Plus  = hamp.Amp0Plus();
  //New terms specific to DM
  fz0Minus = hamp.Ampz0Minus();
  fz0Plus  = hamp.Ampz0Plus();
}
//____________________________________________________________________________
void RSHelicityAmplDM::Print(ostream & stream) const
{
  stream << endl;
  stream << " f(-1) = " << fMinus1 << endl;
  stream << " f(+1) = " << fPlus1  << endl;
  stream << " f(-3) = " << fMinus3 << endl;
  stream << " f(+3) = " << fPlus3  << endl;
  stream << " f(0-) = " << f0Minus << endl;
  stream << " f(0+) = " << f0Plus  << endl;
  //New terms specific to DM
  stream << " fz(0-) = " << fz0Minus << endl;
  stream << " fz(0+) = " << fz0Plus  << endl;
}
//____________________________________________________________________________
void RSHelicityAmplDM::Init(void)
{
  fMinus1 = 0.0;
  fPlus1  = 0.0;
  fMinus3 = 0.0;
  fPlus3  = 0.0;
  f0Minus = 0.0;
  f0Plus  = 0.0;
  //New terms specific to DM
  fz0Minus = 0.0;
  fz0Plus  = 0.0;
}
//____________________________________________________________________________

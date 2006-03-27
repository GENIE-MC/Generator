//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmpl

\brief    A class holding the Rein-Seghal's helicity amplitudes.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          RSHelicityAmplModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "ReinSeghal/RSHelicityAmpl.h"

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



//____________________________________________________________________________
/*!

\class    genie::BenharSpectralFunc

\brief    A realistic spectral function computed using the Local Density
          Aproximation.
          Is a concrete implementation of the SpectralFunctionI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/
//____________________________________________________________________________

#ifndef _BENHAR_SPECTRAL_FUNCTION_H_
#define _BENHAR_SPECTRAL_FUNCTION_H_

#include "Nuclear/NuclearModelI.h"

namespace genie {

class BenharSpectralFunc : public NuclearModelI {

public:
  BenharSpectralFunc();
  BenharSpectralFunc(string config);
  virtual ~BenharSpectralFunc();

  //-- implement the NuclearModelI interface
  bool     GenerateNucleon (const Target & t) const;
  double   Prob            (double p, double w, const Target & t) const;
};

}      // genie namespace
#endif // _BENHAR_SPECTRAL_FUNCTION_H_

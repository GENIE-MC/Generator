//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionI

\brief    Pure abstract base class. Defines the SpectralFunctionI interface
          to be implemented by any algorithmic class implementing a spectral
          function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_I_H_
#define _SPECTRAL_FUNCTION_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class SpectralFunctionI : public Algorithm {

public:

  virtual ~SpectralFunctionI();

  //-- define SpectralFunctionI interface
  
  virtual double Prob(double P_nucleon, double E_nucleus) = 0;

protected:

  SpectralFunctionI();
  SpectralFunctionI(const char * param_set);
};

}      // genie namespace

#endif // _SPECTRAL_FUNCTION_I_H_

//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionFG

\brief    A Fermi Gas equivalent spectral function.

          Is a concrete implementation of the SpectralFunctionI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_FG_H_
#define _SPECTRAL_FUNCTION_FG_H_

#include "Nuclear/SpectralFunctionI.h"

namespace genie {

class SpectralFunctionFG : public SpectralFunctionI {

public:

  SpectralFunctionFG();
  SpectralFunctionFG(string config);  
  ~SpectralFunctionFG();

  //-- SpectralFunctionI interface implementation
  
  double Prob(double P_nucleon, double E_nucleus);
};

}      // genie namespace

#endif // _SPECTRAL_FUNCTION_FG_H_

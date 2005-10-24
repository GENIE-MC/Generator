//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionLDA

\brief    A realistic spectral function computed using the Local Density
          Aproximation.

          Is a concrete implementation of the SpectralFunctionI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_LDA_H_
#define _SPECTRAL_FUNCTION_LDA_H_

#include "Nuclear/SpectralFunctionI.h"

namespace genie {

class SpectralFunctionLDA : public SpectralFunctionI {

public:

  SpectralFunctionLDA();
  SpectralFunctionLDA(string config);
  ~SpectralFunctionLDA();

  //-- SpectralFunctionI interface implementation
  
  double Prob(double P_nucleon, double E_nucleus);
};

}      // genie namespace

#endif // _SPECTRAL_FUNCTION_LDA_H_

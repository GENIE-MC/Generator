//____________________________________________________________________________
/*!

\class    genie::FGMSpectralFunc

\brief    A Fermi Gas equivalent spectral function.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/
//____________________________________________________________________________

#ifndef _FGM_SPECTRAL_FUNC_H_
#define _FGM_SPECTRAL_FUNC_H_

#include "Nuclear/NuclearModelI.h"

namespace genie {

class FGMSpectralFunc : public NuclearModelI {

public:
  FGMSpectralFunc();
  FGMSpectralFunc(string config);  
 ~FGMSpectralFunc();

  //-- implement the NuclearlI interface
  bool     GenerateNucleon (const Target & t) const;
  double   Prob            (double p, double w, const Target & t) const;
};

}      // genie namespace
#endif // _FGM_SPECTRAL_FUNC_H_

//____________________________________________________________________________
/*!

\class    genie::BenharSpectralFunc1D

\brief    
          Implements the NuclearModelI interface.

\ref      

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _BENHAR_SPECTRAL_FUNCTION_1D_H_
#define _BENHAR_SPECTRAL_FUNCTION_1D_H_

#include "Nuclear/NuclearModelI.h"

namespace genie {

class BenharSpectralFunc1D : public NuclearModelI {

public:
  BenharSpectralFunc1D();
  BenharSpectralFunc1D(string config);
  virtual ~BenharSpectralFunc1D();

  //-- implement the NuclearModelI interface
  bool     GenerateNucleon (const Target & t) const;
  double   Prob            (double p, double w, const Target & t) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig();
  TH1D * fHisto;
};

}         // genie namespace
#endif    // _BENHAR_SPECTRAL_FUNCTION_1D_H_


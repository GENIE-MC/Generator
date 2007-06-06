//____________________________________________________________________________
/*!

\class    genie::SpectralFunc1d

\brief    Simpler approach to using spectral functions.
          Implements the NuclearModelI interface.

\ref      

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_1D_H_
#define _SPECTRAL_FUNCTION_1D_H_

#include "Nuclear/NuclearModelI.h"

namespace genie {

class Spline;
class SpectralFunc1d : public NuclearModelI {

public:
  SpectralFunc1d();
  SpectralFunc1d(string config);
  virtual ~SpectralFunc1d();

  //-- implement the NuclearModelI interface
  bool     GenerateNucleon (const Target & t) const;
  double   Prob            (double p, double w, const Target & t) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void     LoadConfig            (void);
  Spline * SelectMomentumDistrib (const Target & target) const;
  double   MaxProb               (const Target & target) const;

  Spline * fSfC12_k;     ///< Benhar C12  spectral func integrated over removal energy
  Spline * fSfFe56_k;    ///< Benhar Fe56 spectral func integrated over removal energy
  double   fMaxC12Prob;  ///<
  double   fMaxFe56Prob; ///<
};

}         // genie namespace
#endif    // _SPECTRAL_FUNCTION_1D_H_


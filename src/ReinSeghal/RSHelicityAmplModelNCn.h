//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelNCn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_NC_N_H_
#define _HELICITY_AMPL_MODEL_NC_N_H_

#include "ReinSeghal/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelNCn : public RSHelicityAmplModelI {

public:

  RSHelicityAmplModelNCn();
  RSHelicityAmplModelNCn(const char * param_set);
  virtual ~RSHelicityAmplModelNCn();

  //-- SPPHelicityAmplModelI interface implementation

  double AmpMinus1 (const Interaction * interaction, const FKR & fkr) const;
  double AmpPlus1  (const Interaction * interaction, const FKR & fkr) const;
  double AmpMinus3 (const Interaction * interaction, const FKR & fkr) const;
  double AmpPlus3  (const Interaction * interaction, const FKR & fkr) const;
  double Amp0Minus (const Interaction * interaction, const FKR & fkr) const;
  double Amp0Plus  (const Interaction * interaction, const FKR & fkr) const;
};

}        // genie namespace

#endif   // _HELICITY_AMPL_MODEL_NC_N_H_

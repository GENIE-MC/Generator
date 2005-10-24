//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the SPPHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_CC_H_
#define _HELICITY_AMPL_MODEL_CC_H_

#include "ReinSeghal/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelCC : public RSHelicityAmplModelI {

public:

  RSHelicityAmplModelCC();
  RSHelicityAmplModelCC(string config);
  virtual ~RSHelicityAmplModelCC();

  //-- SPPHelicityAmplModelI interface implementation

  double AmpMinus1 (const Interaction * interaction, const FKR & fkr) const;
  double AmpPlus1  (const Interaction * interaction, const FKR & fkr) const;
  double AmpMinus3 (const Interaction * interaction, const FKR & fkr) const;
  double AmpPlus3  (const Interaction * interaction, const FKR & fkr) const;
  double Amp0Minus (const Interaction * interaction, const FKR & fkr) const;
  double Amp0Plus  (const Interaction * interaction, const FKR & fkr) const;
};

}        // namespace

#endif   // _HELICITY_AMPL_MODEL_CC_H_



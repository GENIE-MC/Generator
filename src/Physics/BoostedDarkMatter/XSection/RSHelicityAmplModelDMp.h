//____________________________________________________________________________
/*
Dark Matter Resonant Scattering Amplitude Calculation for proton Interactions (header)

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________


#ifndef _HELICITY_AMPL_MODEL_DM_P_H_
#define _HELICITY_AMPL_MODEL_DM_P_H_

#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMI.h"

namespace genie {

class RSHelicityAmplModelDMp : public RSHelicityAmplModelDMI {

public:
  RSHelicityAmplModelDMp();
  RSHelicityAmplModelDMp(string config);
  virtual ~RSHelicityAmplModelDMp();


  const RSHelicityAmplDM & Compute(Resonance_t res, const FKRDM & fkrdm) const;

private:
  mutable RSHelicityAmplDM fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_DM_N_H_

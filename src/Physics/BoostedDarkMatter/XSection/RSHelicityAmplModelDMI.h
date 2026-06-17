//____________________________________________________________________________
/*
RSHelicityAmplModelDMI: the interface for the Dark Matter helicity Amplitudes (header)

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_
#define _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/BoostedDarkMatter/XSection/FKRDM.h"
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplDM.h"

namespace genie {

class RSHelicityAmplModelDMI : public Algorithm
{
public:
//______________________ added for quark charges
  void Configure(string config);
//______________________

  virtual ~RSHelicityAmplModelDMI();

  virtual const RSHelicityAmplDM & Compute(Resonance_t res, const FKRDM & fkrdm) const = 0;



protected:
  RSHelicityAmplModelDMI();
  RSHelicityAmplModelDMI(string name);
  RSHelicityAmplModelDMI(string name, string config);

  //______________added for quark charges
    void LoadConfig(void);

    double fQuV;
    double fQuA;
    double fQdV;
    double fQdA;
  //________________

};

}        // namespace

#endif   // _REIN_SEHGAL_HELICITY_AMPL_MODEL_DM_I_H_

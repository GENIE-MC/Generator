//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelNCp

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free protons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_NC_P_H_
#define _HELICITY_AMPL_MODEL_NC_P_H_

#include "ReinSeghal/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelNCp : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelNCp();
  RSHelicityAmplModelNCp(string config);
  virtual ~RSHelicityAmplModelNCp();

  //-- RSHelicityAmplModelI interface implementation
  RSHelicityAmpl * Compute(Resonance_t res, const FKR & fkr) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  double fSin28w;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_NC_P_H_

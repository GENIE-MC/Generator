//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelNCn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_NC_N_H_
#define _HELICITY_AMPL_MODEL_NC_N_H_

#include "ReinSeghal/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelNCn : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelNCn();
  RSHelicityAmplModelNCn(string config);
  virtual ~RSHelicityAmplModelNCn();

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
#endif   // _HELICITY_AMPL_MODEL_NC_N_H_

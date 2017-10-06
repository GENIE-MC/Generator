//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmplModelNCn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Sehgal's paper.

          Concrete implementation of the MAIDHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_NC_N_H_
#define _MAID_HELICITY_AMPL_MODEL_NC_N_H_

#include "eresMAID/MAIDHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelNCn : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelNCn();
  MAIDHelicityAmplModelNCn(string config);
  virtual ~MAIDHelicityAmplModelNCn();

  // MAIDHelicityAmplModelI interface implementation
  const MAIDHelicityAmpl & Compute(Resonance_t res, const FKR & fkr, const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  mutable MAIDHelicityAmpl fAmpl;
  double fSin28w;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_NC_N_H_

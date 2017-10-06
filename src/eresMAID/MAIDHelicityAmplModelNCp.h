//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmplModelNCp

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free protons, as computed in the Rein-Sehgal's paper.

          Concrete implementation of the MAIDHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_NC_P_H_
#define _MAID_HELICITY_AMPL_MODEL_NC_P_H_

#include "eresMAID/MAIDHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelNCp : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelNCp();
  MAIDHelicityAmplModelNCp(string config);
  virtual ~MAIDHelicityAmplModelNCp();

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
#endif   // _MAID_HELICITY_AMPL_MODEL_NC_P_H_

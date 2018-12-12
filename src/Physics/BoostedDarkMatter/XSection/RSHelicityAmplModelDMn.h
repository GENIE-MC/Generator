//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelDMn

\brief    The Helicity Amplitudes, for all baryon resonances, for NC neutrino
          interactions on free neutrons, as computed in the Rein-Sehgal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_DM_N_H_
#define _HELICITY_AMPL_MODEL_DM_N_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelDMn : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelDMn();
  RSHelicityAmplModelDMn(string config);
  virtual ~RSHelicityAmplModelDMn();

  // RSHelicityAmplModelI interface implementation
  const RSHelicityAmpl & Compute(Resonance_t res, const FKR & fkr) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  mutable RSHelicityAmpl fAmpl;

  double fQVu;
  double fQAu;
  double fQVd;
  double fQAd;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_DM_N_H_

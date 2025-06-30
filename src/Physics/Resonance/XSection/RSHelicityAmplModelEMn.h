//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelEMn

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free neutrons, as computed in the
          Rein-Sehgal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  March 30, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_EM_N_H_
#define _HELICITY_AMPL_MODEL_EM_N_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelEMn : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelEMn();
  RSHelicityAmplModelEMn(string config);
  virtual ~RSHelicityAmplModelEMn();

  // RSHelicityAmplModelI interface implementation
  const RSHelicityAmpl & Compute(Resonance_t res, const FKR & fkr) const;

private:
  mutable RSHelicityAmpl fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_EM_N_H_

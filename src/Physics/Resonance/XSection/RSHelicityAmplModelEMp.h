//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelEMp

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free protons, as computed in the
          Rein-Sehgal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  March 30, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_EM_P_H_
#define _HELICITY_AMPL_MODEL_EM_P_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelEMp : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelEMp();
  RSHelicityAmplModelEMp(string config);
  virtual ~RSHelicityAmplModelEMp();

  // RSHelicityAmplModelI interface implementation
  const RSHelicityAmpl & Compute(Resonance_t res, const FKR & fkr) const;

private:
  mutable RSHelicityAmpl fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_EM_N_H_

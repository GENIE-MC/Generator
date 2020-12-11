//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Sehgal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_CC_H_
#define _HELICITY_AMPL_MODEL_CC_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

namespace genie {

class RSHelicityAmplModelCC : public RSHelicityAmplModelI {

public:
  RSHelicityAmplModelCC();
  RSHelicityAmplModelCC(string config);
  virtual ~RSHelicityAmplModelCC();

  // RSHelicityAmplModelI interface implementation
 const RSHelicityAmpl & Compute(Resonance_t res, const FKR & fkr) const;

private:
  mutable RSHelicityAmpl fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_CC_H_

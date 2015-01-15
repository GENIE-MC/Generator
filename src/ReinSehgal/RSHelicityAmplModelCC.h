//____________________________________________________________________________
/*!

\class    genie::RSHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Seghal's paper.

          Concrete implementation of the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_CC_H_
#define _HELICITY_AMPL_MODEL_CC_H_

#include "ReinSeghal/RSHelicityAmplModelI.h"

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



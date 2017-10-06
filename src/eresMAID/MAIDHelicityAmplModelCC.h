//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmplModelCC

\brief    The Helicity Amplitudes, for all baryon resonances, for CC neutrino
          interactions on free nucleons, as computed in the Rein-Sehgal's paper.

          Concrete implementation of the MAIDHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_CC_H_
#define _MAID_HELICITY_AMPL_MODEL_CC_H_

#include "eresMAID/MAIDHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelCC : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelCC();
  MAIDHelicityAmplModelCC(string config);
  virtual ~MAIDHelicityAmplModelCC();

  // MAIDHelicityAmplModelI interface implementation
 const MAIDHelicityAmpl & Compute(Resonance_t res, const FKR & fkr, const Interaction * i) const;

private:
  mutable MAIDHelicityAmpl fAmpl; 
};

}        // genie namespace
#endif   // _MAID_HELICITY_AMPL_MODEL_CC_H_



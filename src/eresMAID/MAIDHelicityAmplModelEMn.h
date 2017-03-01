//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmplModelEMn

\brief    The Helicity Amplitudes, for all baryon resonances, for Electro-
          Magnetic (EM) interactions on free neutrons, as computed in the
          Rein-Sehgal's paper.

          Concrete implementation of the MAIDHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 30, 2005

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HELICITY_AMPL_MODEL_EM_N_H_
#define _HELICITY_AMPL_MODEL_EM_N_H_

#include "eresMAID/MAIDHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelEMn : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelEMn();
  MAIDHelicityAmplModelEMn(string config);
  virtual ~MAIDHelicityAmplModelEMn();

  // MAIDHelicityAmplModelI interface implementation
  const MAIDHelicityAmpl & Compute(Resonance_t res, const FKR & fkr, const Interaction * i) const;

private:
  mutable MAIDHelicityAmpl fAmpl;
};

}        // genie namespace
#endif   // _HELICITY_AMPL_MODEL_EM_N_H_

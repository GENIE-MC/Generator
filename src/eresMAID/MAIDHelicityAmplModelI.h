//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmplModelI

\brief    Pure abstract base class. Defines the MAIDHelicityAmplModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  July 10, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REIN_SEHGAL_HELICITY_AMPL_MODEL_I_H_
#define _REIN_SEHGAL_HELICITY_AMPL_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "BaryonResonance/BaryonResonance.h"
#include "eresMAID/FKR.h"
#include "eresMAID/MAIDHelicityAmpl.h"

namespace genie {

class MAIDHelicityAmplModelI : public Algorithm
{
public:
  virtual ~MAIDHelicityAmplModelI();

  // define the MAIDHelicityAmplModelI interface
  virtual const MAIDHelicityAmpl & Compute(Resonance_t res, const FKR & fkr, const Interaction * interaction) const = 0;

protected:
  MAIDHelicityAmplModelI();
  MAIDHelicityAmplModelI(string name);
  MAIDHelicityAmplModelI(string name, string config);
};

}        // namespace

#endif   // _REIN_SEHGAL_HELICITY_AMPL_MODEL_I_H_




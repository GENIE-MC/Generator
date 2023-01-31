//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmplModelI

\brief    Pure abstract base class. Defines the MAIDHelicityAmplModelI interface.

\author   Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
          Tel Aviv University

\created  January 2023

\cpright  Copyright (c) 2023-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_I_H_
#define _MAID_HELICITY_AMPL_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"

namespace genie {

class MAIDHelicityAmplModelI : public Algorithm
{
public:
  virtual ~MAIDHelicityAmplModelI();

  // define the MAIDHelicityAmplModelI interface
  virtual const MAIDHelicityAmpl & Compute( const Interaction * interaction ) const = 0;

protected:
  MAIDHelicityAmplModelI();
  MAIDHelicityAmplModelI(string name);
  MAIDHelicityAmplModelI(string name, string config);
};

}        // namespace

#endif   // _MAID_HELICITY_AMPL_MODEL_I_H_

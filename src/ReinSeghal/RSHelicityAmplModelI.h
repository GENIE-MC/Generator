//____________________________________________________________________________
/*!
\class    genie::RSHelicityAmplModelI

\brief    Pure abstract base class. Defines the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 10, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_
#define _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "BaryonResonance/BaryonResonance.h"
#include "ReinSeghal/FKR.h"

namespace genie {

class RSHelicityAmpl;
class RSHelicityAmplModelI : public Algorithm
{
public:
  virtual ~RSHelicityAmplModelI();

  //-- define the RSHelicityAmplModelI interface
  virtual RSHelicityAmpl * Compute(Resonance_t res, const FKR & fkr) const = 0;

protected:
  RSHelicityAmplModelI();
  RSHelicityAmplModelI(string name);
  RSHelicityAmplModelI(string name, string config);
};

}        // namespace

#endif   // _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_




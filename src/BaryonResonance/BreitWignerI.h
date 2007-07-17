//____________________________________________________________________________
/*!

\class    genie::BreitWignerI

\brief    Pure abstract base class. Defines the BreitWignerI interface to
          be implemented by any algorithmic class modeling a Breit Wigner
          function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________


#ifndef _BREIT_WIGNER_I_H_
#define _BREIT_WIGNER_I_H_

#include "Algorithm/Algorithm.h"
#include "BaryonResonance/BaryonResonance.h"

namespace genie {

class BreitWignerI : public Algorithm {

public:
  virtual ~BreitWignerI();

  //! Evaluate the Breit-Wigner function for the input resonance at the input mass
  virtual double Eval(Resonance_t res, double W) const = 0;

protected:
  BreitWignerI();
  BreitWignerI(string name);
  BreitWignerI(string name, string config);
};

}         // genie namespace

#endif    // _BREIT_WIGNER_I_H_

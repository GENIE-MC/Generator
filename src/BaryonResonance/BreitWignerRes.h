//____________________________________________________________________________
/*!

\class    genie::BreitWignerRes

\brief    Concrete implementation of the BreitWignerI interface:
          Simple Breit-Wigner distribution with no L-dependent width.

          It is similar with the breit_wigner function but rather than
          specifying the Breit-Wigner parameters directly, you specify a
          resonance name and the concrete implementation of BaryonResDataSetI
          to be looked up for extracting those parameters.

          Pre-configured instances can be obtained from the AlgFactory.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_RES_H_
#define _BREIT_WIGNER_RES_H_

#include "BaryonResonance/BreitWignerI.h"

namespace genie {

class BreitWignerRes : public BreitWignerI {

public:
  BreitWignerRes(); 
  BreitWignerRes(const char * param_set); 
  virtual ~BreitWignerRes(); 

  double Eval(double W) const;
};

}         // genie namespace

#endif    // _BREIT_WIGNER_RES_H_

//____________________________________________________________________________
/*!

\class    genie::BreitWignerLRes

\brief    Concrete implementation of the BreitWignerI interface:
          A realistic Breit-Wigner distribution with L-dependent width.
          
          It is similar with the breit_wigner_L function but rather than 
          specifying the Breit-Wigner parameters directly, you specify a 
          resonance name and the concrete implementation of BaryonResDataSetI 
          to be looked up for extracting those parameters.

          Pre-configured instances can be obtained from the AlgFactory
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_L_RES_H_
#define _BREIT_WIGNER_L_RES_H_

#include "BaryonResonance/BreitWignerI.h"

namespace genie {

class BreitWignerLRes : public BreitWignerI {

public:

  BreitWignerLRes();
  BreitWignerLRes(string config);
  ~BreitWignerLRes(); 

  double Eval(Double_t W) const;
};

}        // genie namespace

#endif   // _BREIT_WIGNER_L_RES_H_

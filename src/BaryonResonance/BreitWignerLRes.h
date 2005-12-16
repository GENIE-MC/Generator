//____________________________________________________________________________
/*!

\class    genie::BreitWignerLRes

\brief    Concrete implementation of the BreitWignerI interface:
          A realistic Breit-Wigner distribution with L-dependent width.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_L_RES_H_
#define _BREIT_WIGNER_L_RES_H_

#include "BaryonResonance/BreitWignerI.h"

namespace genie {

class BaryonResDataSetI;

class BreitWignerLRes : public BreitWignerI {

public:

  BreitWignerLRes();
  BreitWignerLRes(string config);
  ~BreitWignerLRes();

  //-- implement the BreitWignerI interface
  double Eval(Resonance_t res, double W) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  //-- load sub-algorithm specified by configuration option
  void LoadSubAlg(void);

  const BaryonResDataSetI * fBaryonResDataSet;
};

}        // genie namespace

#endif   // _BREIT_WIGNER_L_RES_H_

//____________________________________________________________________________
/*!

\class    genie::BreitWignerRes

\brief    Concrete implementation of the BreitWignerI interface:
          Simple Breit-Wigner distribution with no L-dependent width.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 22, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_RES_H_
#define _BREIT_WIGNER_RES_H_

#include "BaryonResonance/BreitWignerI.h"

namespace genie {

class BaryonResDataSetI;

class BreitWignerRes : public BreitWignerI {

public:
  BreitWignerRes();
  BreitWignerRes(string config);
  virtual ~BreitWignerRes();

  //-- implement the BreitWignerI interface
  double Eval(Resonance_t res, double W) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  //-- load sub-algorithm specified by configuration option
  void LoadConfig(void);

  //-- private data members
  const BaryonResDataSetI * fBaryonResDataSet;
};

}         // genie namespace

#endif    // _BREIT_WIGNER_RES_H_


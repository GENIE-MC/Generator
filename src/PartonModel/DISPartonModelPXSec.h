//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    DIS differential (d^2xsec/dxdy) cross section

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _DIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Base/DISStructureFunc.h"

namespace genie {

class DISStructureFuncModelI;
class MultiplicityProbModelI;

class DISPartonModelPXSec : public XSecAlgorithmI {

public:
  DISPartonModelPXSec();
  DISPartonModelPXSec(string config);
  virtual ~DISPartonModelPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig                  (void);
  double DISRESJoinSuppressionFactor (const Interaction * in) const;

  mutable DISStructureFunc fDISSF;

  //! configuration data

  bool   fUsingDisResJoin;  ///< use a DIS/RES joining scheme?
  double fWcut;             ///< apply DIS/RES joining scheme < Wcut

  const DISStructureFuncModelI * fDISSFModel;    ///< SF model
  const MultiplicityProbModelI * fMultProbModel; ///< hadronic multip. model
};

}       // genie namespace
#endif  // _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

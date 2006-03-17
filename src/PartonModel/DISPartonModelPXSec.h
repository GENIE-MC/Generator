//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    DIS differential (d^2xsec/dxdy) cross section

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class DISStructureFuncModelI;

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
  void LoadSubAlg (void);
  const DISStructureFuncModelI * fDISSFModel;
};

}       // genie namespace
#endif  // _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    Computes DIS differential cross sections.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
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
class HadronizationModelI;
class XSecIntegratorI;

class DISPartonModelPXSec : public XSecAlgorithmI {

public:
  DISPartonModelPXSec();
  DISPartonModelPXSec(string config);
  virtual ~DISPartonModelPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig                  (void);
  double DISRESJoinSuppressionFactor (const Interaction * in) const;

  mutable DISStructureFunc fDISSF;
  bool                     fInInitPhase;

  const DISStructureFuncModelI * fDISSFModel;         ///< SF model
  const HadronizationModelI *    fHadronizationModel; ///< hadronic multip. model
  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  const XSecAlgorithmI * fCharmProdModel;

  bool   fUsingDisResJoin;  ///< use a DIS/RES joining scheme?
  bool   fUseCache;         ///< cache reduction factors used in joining scheme
  double fWcut;             ///< apply DIS/RES joining scheme < Wcut
  double fScale;            ///< cross section scaling factor
};

}       // genie namespace
#endif  // _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

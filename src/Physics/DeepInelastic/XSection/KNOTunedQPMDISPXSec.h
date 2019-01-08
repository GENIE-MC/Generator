//____________________________________________________________________________
/*!

\class    genie::KNOTunedQPMDISPXSec

\brief    Computes DIS differential cross sections.
          Is a concrete implementation of the XSecAlgorithmI interface.


\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KNO_DIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _KNO_DIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DeepInelastic/XSection/QPMDISPXSec.h"
#include "Physics/Hadronization/KNOHadronization.h"


namespace genie {

class HadronizationModelI;

class KNOTunedQPMDISPXSec : public XSecAlgorithmI {

public:
  KNOTunedQPMDISPXSec();
  KNOTunedQPMDISPXSec(string config);
  virtual ~KNOTunedQPMDISPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig                  (void);
  double DISRESJoinSuppressionFactor (const Interaction * in) const;

  const KNOHadronization *    fHadronizationModel; ///< hadronic multip. model
  const QPMDISPXSec *         fDISModel ;
  const XSecIntegratorI *     fXSecIntegrator;     ///< diff. xsec integrator

  bool   fUseCache;         ///< cache reduction factors used in joining scheme
  double fWcut;             ///< apply DIS/RES joining scheme < Wcut
};

}       // genie namespace
#endif  // _KNO_DIS_PARTON_MODEL_PARTIAL_XSEC_H_

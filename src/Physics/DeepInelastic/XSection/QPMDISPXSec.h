//____________________________________________________________________________
/*!

\class    genie::QPMDISPXSec

\brief    Computes DIS differential cross sections.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 05, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _DIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"

namespace genie {

class DISStructureFuncModelI;
class HadronizationModelI;
class XSecIntegratorI;

class QPMDISPXSec : public XSecAlgorithmI {

public:
  QPMDISPXSec();
  QPMDISPXSec(string config);
  virtual ~QPMDISPXSec();

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

  mutable DISStructureFunc fDISSF;
  bool                     fInInitPhase;

  const DISStructureFuncModelI * fDISSFModel;         ///< SF model
  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  const XSecAlgorithmI * fCharmProdModel;

  double fCCScale;            ///< cross section scaling factor
  double fNCScale;            ///< cross section scaling factor
  double fEMScale;            ///< cross section scaling factor
  double fSin48w;           ///< sin^4(Weingberg angle)
};

}       // genie namespace
#endif  // _DIS_PARTON_MODEL_PARTIAL_XSEC_H_

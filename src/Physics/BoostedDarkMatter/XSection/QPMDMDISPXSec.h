//____________________________________________________________________________
/*!

\class    genie::QPMDMDISPXSec

\brief    Computes DMDIS differential cross sections.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.03300

\author   Joshua Berger <jberger \at physics.wisc.edu
          University of Wisconsin-Madison

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 4, 2017

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DMDIS_PARTON_MODEL_PARTIAL_XSEC_H_
#define _DMDIS_PARTON_MODEL_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"

namespace genie {

class DISStructureFuncModelI;
class HadronizationModelI;
class XSecIntegratorI;

class QPMDMDISPXSec : public XSecAlgorithmI {

public:
  QPMDMDISPXSec();
  QPMDMDISPXSec(string config);
  virtual ~QPMDMDISPXSec();

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
  double DMDISRESJoinSuppressionFactor (const Interaction * in) const;

  mutable DISStructureFunc fDISSF;
  bool                     fInInitPhase;

  const DISStructureFuncModelI * fDISSFModel;         ///< SF model
  const HadronizationModelI *    fHadronizationModel; ///< hadronic multip. model
  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  const XSecAlgorithmI * fCharmProdModel;

  bool   fUsingDisResJoin;  ///< use a DMDIS/RES joining scheme?
  bool   fUseCache;         ///< cache reduction factors used in joining scheme
  double fWcut;             ///< apply DMDIS/RES joining scheme < Wcut
  double fScale;            ///< cross section scaling factor
  double fSin48w;           ///< sin^4(Weingberg angle)
  int    fVelMode;          ///< velcoity dependence for xsec
  double fMedMass;          ///< Mediator mass
  double fgzp;              ///< Coupling to the mediator Zprime
};

}       // genie namespace
#endif  // _DMDIS_PARTON_MODEL_PARTIAL_XSEC_H_

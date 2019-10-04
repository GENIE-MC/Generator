//____________________________________________________________________________
/*!

\class    genie::EmpiricalMECPXSec2015

\brief    Computes the MEC differential cross section.
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Steve Dytman <dytman+ \at pitt.edu>
          Pittsburgh University

\created  May 05, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_PXSEC_H_
#define _MEC_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class EmpiricalMECPXSec2015 : public XSecAlgorithmI {

public:
  EmpiricalMECPXSec2015();
  EmpiricalMECPXSec2015(string config);
  virtual ~EmpiricalMECPXSec2015();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig (void);

  double fMq2d;       ///< toy model param: `mass' in dipole (Q2 - dependence) form factor (GeV)
  double fMass;       ///< toy model param: peak  of W distribution (GeV)
  double fWidth;      ///< toy model param: width of W distribution (GeV)
  double fMECAPower;      ///< power of A relative to carbon

  double fFracPN_NC;     ///< toy model param: fraction of nucleon pairs that are pn, not nn or pp
  double fFracPN_CC;     ///< toy model param: fraction of nucleon pairs that are pn, not nn or pp
  double fFracPN_EM;     ///< toy model param: fraction of nucleon pairs that are pn, not nn or pp

  double fFracCCQE;   ///< empirical model param: MEC cross section is taken to be this fraction of CCQE cross section
  double fFracNCQE;   ///< empirical model param: MEC cross section is taken to be this fraction of NCQE cross section
  double fFracEMQE;   ///< empirical model param: MEC cross section is taken to be this fraction of Rosenbluth xs

  const XSecAlgorithmI * fXSecAlgCCQE; ///< cross section algorithm for CCQE
  const XSecAlgorithmI * fXSecAlgNCQE; ///< cross section algorithm for NCQE
  const XSecAlgorithmI * fXSecAlgEMQE; ///< cross section algorithm for EMQE

};

}       // genie namespace
#endif  // _MEC_PARTIAL_XSEC_H_

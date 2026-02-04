//____________________________________________________________________________
/*!

\class    genie::SmithMonizQELCCPXSec

\brief    Computes neutrino-nucleon(nucleus) QELCC differential cross section.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      [1] R.A.Smith and E.J.Moniz,
              Nuclear Physics  B43, (1972) 605-622 \n
          [2] K.S. Kuzmin, V.V. Lyubushkin, V.A.Naumov,
              Eur. Phys. J. C54, (2008) 517-538

\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research \n

          adapted from  fortran code provided by: \n

          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>
          Joint Institute for Nuclear Research \n

          Vladimir Lyubushkin
          Joint Institute for Nuclear Research \n

          Vadim Naumov <vnaumov@theor.jinr.ru>
          Joint Institute for Nuclear Research  \n

          based on code of: \n
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 05, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SMITH_MONITZ_QELCC_CROSS_SECTION_H_
#define _SMITH_MONITZ_QELCC_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"


namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class SmithMonizQELCCPXSec : public XSecAlgorithmI {

public:
  SmithMonizQELCCPXSec();
  SmithMonizQELCCPXSec(string config);
  virtual ~SmithMonizQELCCPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t kps) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  mutable SmithMonizUtils * sm_utils;

  void   LoadConfig (void);
  double d3sQES_dQ2dvdkF_SM (const Interaction * interaction) const;
  double dsQES_dQ2_SM(const Interaction * interaction) const;
  double d2sQES_dQ2dv_SM(const Interaction * i) const;

  double                       fXSecScale;        ///< external xsec scaling factor
  mutable QELFormFactors       fFormFactors;
  const QELFormFactorsModelI * fFormFactorsModel;
  const XSecIntegratorI *      fXSecIntegrator;
  double                       fVud2;             ///< |Vud|^2(square of magnitude ud-element of CKM-matrix)


};


} // genie namespace

#endif  //_SMITH_MONITZ_QELCC_CROSS_SECTION_H_

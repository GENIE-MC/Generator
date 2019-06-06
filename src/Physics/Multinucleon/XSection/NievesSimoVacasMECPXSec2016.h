//____________________________________________________________________________
/*!

\class    genie::NievesSimoVacasMECPXSec2016

\brief    Computes the Valencia MEC model differential cross section.
          Uses precomputed hadon tensor tables.
          Is a concrete implementation of the XSecAlgorithmI interface. 

\author   Code contributed by J. Schwehr, D. Cherdack, R. Gran and described
          in arXiv:1601.02038 and some of the refereces there-in,
          in particular PRD 88 (2013) 113007

          Substantial code refactorizations by the core GENIE group.

\ref      J. Nieves, I. Ruiz Simo, M.J. Vicente Vacas,
          Inclusive quasi-elastic neutrino reactions, PRC 83 (2011) 045501

\created  Mar 22, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NIEVES_SIMO_VACAS_MEC_PXSEC_2016_H_
#define _NIEVES_SIMO_VACAS_MEC_PXSEC_2016_H_

#include <vector>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"

using std::vector;

namespace genie {

class XSecIntegratorI;

class NievesSimoVacasMECPXSec2016 : public XSecAlgorithmI {

public:
  NievesSimoVacasMECPXSec2016();
  NievesSimoVacasMECPXSec2016(string config);
  virtual ~NievesSimoVacasMECPXSec2016();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:

  // Load algorithm configuration
  void LoadConfig (void);

  double                   fXSecScale;        ///< external xsec scaling factor

  const XSecIntegratorI *  fXSecIntegrator; // Numerical integrator (GSL)

//double fQ3Max;
};
  
}       // genie namespace
#endif  // _NIEVES_SIMO_VACAS_MEC_PXSEC_2016_H_

//____________________________________________________________________________
/*!
  \class    genie::BertuzzoDNuCOHPXSec

  \brief    Differential cross section for v+As coherent elastic scattering.Coherent DNu. \n
            Is a concrete implementation of the XSecAlgorithmI interface. \n

  \ref      E.Bertuzzo, S.Jana, P.A.N.Machado, R. Zukanovich Funcal
            PhysRevLett.121.241801 (2018)

  \author   Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
            University of Sussex

            Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

  \created  June 12, 2020

  \cpright  Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _BERTUZZO_DNu_COH_CROSS_SECTION_H_
#define _BERTUZZO_DNu_COH_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/DarkNeutrino/XSection/EngelFormFactor.h"

namespace genie {

class XSecIntegratorI;

class BertuzzoDNuCOHPXSec : public XSecAlgorithmI {

public:
  BertuzzoDNuCOHPXSec();
  BertuzzoDNuCOHPXSec(string config);
  virtual ~BertuzzoDNuCOHPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const override;
  double Integral        (const Interaction * i) const override;
  bool   ValidProcess    (const Interaction * i) const override;
  bool   ValidKinematics (const Interaction * i) const override;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config) override;
  void Configure (string param_set) override;

private:

  void LoadConfig(void);

  const XSecIntegratorI * fXSecIntegrator;  ///< cross section integrator
  const EngelFormFactor * fFF; ///< Engel Form Factor algorithm

  double fEps2;

  std::array<double, 4> fMixing2s;
  // mixing angles between the neutrinos (including a 4rth sterile one)
  // to the dark neutrino

  double fAlpha_D;

  double fDNuMass, fDNuMass2;
  double fDMediatorMass, fDMediatorMass2;

};

} // genie namespace


#endif // _BERTUZZO_DNu_COH_CROSS_SECTION_H_

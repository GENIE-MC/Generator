//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <TMath.h>
#include <iostream>
#include "Conventions/Constants.h"

#include "Algorithm/AlgConfigPool.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "AlvarezRuso/ARConstants.h"

namespace genie {
namespace alvarezruso {
  
  
ARConstants::ARConstants() {
  reg = genie::AlgConfigPool::Instance()->GlobalParameterList();
  COHAR_Ma_Nuc      = reg->GetDouble("COHAR-Ma-Nuc");
  COHAR_Mv_Nuc      = reg->GetDouble("COHAR-Mv-Nuc");
  COHAR_Ma_Delta    = reg->GetDouble("COHAR-Ma-Delta");
  COHAR_Mv_Delta    = reg->GetDouble("COHAR-Mv-Delta");
  COHAR_GA0         = reg->GetDouble("COHAR-GA0");
  COHAR_Rho0        = reg->GetDouble("COHAR-Rho0");
  COHAR_a4          = reg->GetDouble("COHAR-a4");
  COHAR_a5          = reg->GetDouble("COHAR-a5");
  COHAR_b4          = reg->GetDouble("COHAR-b4");
  COHAR_b5          = reg->GetDouble("COHAR-b5");
  COHAR_fPi_byHbar  = reg->GetDouble("COHAR-fPi") / HBar();
  COHAR_fStar       = reg->GetDouble("COHAR-fStar");
  fCosCabibboAngle  = TMath::Cos(reg->GetDouble("CabibboAngle"));
  fSinWeinbergAngle = TMath::Sin(reg->GetDouble("WeinbergAngle"));
  
  massElectron = genie::PDGLibrary::Instance()->Find(genie::kPdgElectron)->Mass() / HBar();
  massMuon     = genie::PDGLibrary::Instance()->Find(genie::kPdgMuon)->Mass() / HBar();
  massTau      = genie::PDGLibrary::Instance()->Find(genie::kPdgTau)->Mass() / HBar();
  massProton   = genie::PDGLibrary::Instance()->Find(genie::kPdgProton)->Mass() / HBar();
  massNeutron  = genie::PDGLibrary::Instance()->Find(genie::kPdgNeutron)->Mass() / HBar();
  massNucleon  = (massProton + massNeutron)/2.0;
  massNucleon2 = massNucleon*massNucleon;
  massDeltaP   = genie::PDGLibrary::Instance()->Find(genie::kPdgP33m1232_DeltaP)->Mass() / HBar();
  massDelta0   = genie::PDGLibrary::Instance()->Find(genie::kPdgP33m1232_Delta0)->Mass() / HBar();
  massPiP      = genie::PDGLibrary::Instance()->Find(genie::kPdgPiP)->Mass() / HBar();
  massPi0      = genie::PDGLibrary::Instance()->Find(genie::kPdgPi0)->Mass() / HBar();
  
  ncFactor = 1.0 - 2.0*fSinWeinbergAngle*fSinWeinbergAngle;
}

ARConstants::~ARConstants() {
}

double ARConstants::HBar() {
  // Alvarez-Ruso model is expressed in units of fermi. Need a
  // conversion from GeV to fm for energy variables.
  return 0.19733;
}
double ARConstants::Ma_Nucleon() {
  return COHAR_Ma_Nuc;
}
double ARConstants::Mv_Nucleon() {
  return COHAR_Mv_Nuc;
}
double ARConstants::Ma_Delta() {
  return COHAR_Ma_Delta;
}
double ARConstants::Mv_Delta() {
  return COHAR_Mv_Delta;
}
double ARConstants::GAxial() {
  return COHAR_GA0;
}
double ARConstants::Rho0() {
  return COHAR_Rho0;
}
double ARConstants::CA4_A() {
  return COHAR_a4;
}
double ARConstants::CA5_A() {
  return COHAR_a5;
}
double ARConstants::CA4_B() {
  return COHAR_b4;
}
double ARConstants::CA5_B() {
  return COHAR_b5;
}
double ARConstants::PiDecayConst() {
  return COHAR_fPi_byHbar;
}
double ARConstants::DeltaNCoupling() {
  return COHAR_fStar;
}
double ARConstants::CosCabibboAngle() {
  return fCosCabibboAngle;
}
double ARConstants::SinWeinbergAngle() {
  return fSinWeinbergAngle;
}
double ARConstants::GFermi() {
  return (genie::constants::kGF * HBar() * HBar());
}
double ARConstants::ElectronMass() {
  return massElectron;
}
double ARConstants::MuonMass() {
  return massMuon;
}
double ARConstants::TauMass() {
  return massTau;
}
double ARConstants::ProtonMass() {
  return massProton;
}
double ARConstants::NeutronMass() {
  return massNeutron;
}
double ARConstants::NucleonMass() {
  return massNucleon;
}
double ARConstants::NucleonMassSq() {
  return massNucleon2;
}
double ARConstants::DeltaPMass() {
  return massDeltaP;
}
double ARConstants::Delta0Mass() {
  return massDelta0;
}
double ARConstants::PiPMass() {
  return massPiP;
}
double ARConstants::Pi0Mass() {
  return massPi0;
}
double ARConstants::cm38Conversion() {
  return 1E12;
}

double ARConstants::NCFactor() {
  return ncFactor;
}

} //namespace alvarezruso

} //namespace genie

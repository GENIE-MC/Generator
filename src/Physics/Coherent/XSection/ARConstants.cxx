//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include <TMath.h>
#include <iostream>
#include "Framework/Conventions/Constants.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/Coherent/XSection/ARConstants.h"

namespace genie {
namespace alvarezruso {
  
  
ARConstants::ARConstants() {


  COHAR_Ma_Nuc      =  1.000 ;
  COHAR_Mv_Nuc      =  0.840 ;
  COHAR_Ma_Delta    =  1.280 ;
  COHAR_Mv_Delta    =  0.730 ;
  COHAR_GA0         =  1.2670 ;
  COHAR_Rho0        =  0.17 ;
  COHAR_a4          =  -1.21 ;
  COHAR_a5          =  -1.21 ;
  COHAR_b4          = 2.0 ;
  COHAR_b5          = 2.0 ;
  COHAR_fPi_byHbar  = 0.093 / HBar();
  COHAR_fStar       = 2.13 ;
  fCosCabibboAngle  = TMath::Cos( 0.22853207 ) ;
  fSinWeinbergAngle = TMath::Sin( 0.49744211 ) ;
  
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

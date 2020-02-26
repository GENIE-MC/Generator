//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TMath.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithFFNC.h"

using namespace genie;
using namespace genie::constants;

namespace {

  // Helper function that sets a sign based on the isospin of the hit nucleon.
  // This function chooses + (-) for proton (neutron).
  int nucleon_sign( const Interaction* interaction ) {
    int sign = 1;
    int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
    if ( pdg::IsNeutron(hit_nuc_pdg) ) sign = -1;
    return sign;
  }

}

//____________________________________________________________________________
LwlynSmithFFNC::LwlynSmithFFNC() :
LwlynSmithFF("genie::LwlynSmithFFNC")
{

}
//____________________________________________________________________________
LwlynSmithFFNC::LwlynSmithFFNC(string config) :
LwlynSmithFF("genie::LwlynSmithFFNC", config)
{

}
//____________________________________________________________________________
LwlynSmithFFNC::~LwlynSmithFFNC()
{

}
//____________________________________________________________________________
void LwlynSmithFFNC::LoadConfig(void)
{
  // First configure the class members inherited from LwlynSmithFF normally
  LwlynSmithFF::LoadConfig();

  // Then configure the NC strange form factors model
  fStrangeNCFFModel = dynamic_cast< const NCELStrangeFormFactorsModelI* >(
    this->SubAlg("StrangeFormFactorsModel") );
  assert( fStrangeNCFFModel );
}
//____________________________________________________________________________
double LwlynSmithFFNC::F1V(const Interaction* interaction) const
{
  int sign = nucleon_sign( interaction );

  // Calculate F1V-CC
  double F1V_CC = LwlynSmithFF::F1V( interaction );

  // Calculate F2V-EM
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  double F1V_EM = 0.;
  if ( pdg::IsNeutron(hit_nuc_pdg) ) F1V_EM = this->F1N( interaction );
  else F1V_EM = this->F1P( interaction );

  // Calculate the strange NC form factor
  double F1Vs = fStrangeNCFFModel->F1Vs( interaction );

  double F1V_NC = sign*0.5*F1V_CC - 2.*fSin28w*F1V_EM - 0.5*F1Vs;
  return F1V_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::xiF2V(const Interaction* interaction) const
{
  int sign = nucleon_sign( interaction );

  // Calculate F2V_CC
  double F2V_CC = LwlynSmithFF::xiF2V( interaction );

  // Calculate F2V_EM
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  double F2V_EM = 0.;
  if ( pdg::IsNeutron(hit_nuc_pdg) ) F2V_EM = this->F2N( interaction );
  else F2V_EM = this->F2P( interaction );

  // Calculate the strange NC form factor
  double F2Vs = fStrangeNCFFModel->F2Vs( interaction );

  // Calculate F2V-NC
  double F2V_NC = sign*0.5*F2V_CC - 2.*fSin28w*F2V_EM - 0.5*F2Vs;
  return F2V_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::FA(const Interaction* interaction) const
{
  int sign = nucleon_sign( interaction );

  // Calculate FA_CC
  double FA_CC = LwlynSmithFF::FA( interaction );

  // Calculate the strange NC form factor
  double FAs = fStrangeNCFFModel->FAs( interaction );

  // Calculate FA_NC
  double FA_NC = sign*0.5*FA_CC - 0.5*FAs;
  return FA_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::Fp(const Interaction * interaction) const
{
  int sign = nucleon_sign( interaction );

  // Calculate FP_CC
  double FP_CC = LwlynSmithFF::Fp( interaction );

  // Calculate the strange NC form factor
  double FPs = fStrangeNCFFModel->FPs( interaction );

  // Calculate FP_NC
  double FP_NC = sign*0.5*FP_CC - 0.5*FPs;
  return FP_NC;
}
//____________________________________________________________________________

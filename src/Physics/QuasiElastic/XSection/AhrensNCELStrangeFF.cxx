 //____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/XSection/AhrensNCELStrangeFF.h"

using namespace genie;

//____________________________________________________________________________
AhrensNCELStrangeFF::AhrensNCELStrangeFF()
  : NCELStrangeFormFactorsModelI("genie::AhrensNCELStrangeFF")
{

}
//____________________________________________________________________________
AhrensNCELStrangeFF::AhrensNCELStrangeFF(string config)
  : NCELStrangeFormFactorsModelI("genie::AhrensNCELStrangeFF", config)
{

}
//____________________________________________________________________________
AhrensNCELStrangeFF::~AhrensNCELStrangeFF()
{

}
//____________________________________________________________________________
void AhrensNCELStrangeFF::Configure(const Registry& config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensNCELStrangeFF::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void AhrensNCELStrangeFF::LoadConfig(void)
{
  fAxFFModel = dynamic_cast< const AxialFormFactorModelI* >(
    this->SubAlg("AxialFormFactorModel") );

  assert( fAxFFModel );
  fAxFF.SetModel( fAxFFModel );

  // TODO: document this (cite Ahrens paper)
  GetParam( "EL-Axial-Eta", fEta ) ;
}

double AhrensNCELStrangeFF::FAs(const Interaction* interaction) const
{
  // Compute FA-CC
  fAxFF.Calculate( interaction );
  double FA_CC = fAxFF.FA();

  int sign = 1;
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  if ( pdg::IsNeutron(hit_nuc_pdg) ) sign = -1;

  // Compute the strange axial NC form factor
  double FAs = -1. * fEta * FA_CC;
  return FAs;
}

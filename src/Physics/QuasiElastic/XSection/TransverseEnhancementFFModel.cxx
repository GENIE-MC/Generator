//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Brian Coopersmith, University of Rochester

 For the class documentation see the corresponding header file.              

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/TransverseEnhancementFFModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Utils/ConfigIsotopeMapUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::config;

//____________________________________________________________________________
TransverseEnhancementFFModel::TransverseEnhancementFFModel() :
ELFormFactorsModelI("genie::TransverseEnhancementFFModel")
{

}
//____________________________________________________________________________
TransverseEnhancementFFModel::TransverseEnhancementFFModel(string config) :
ELFormFactorsModelI("genie::TransverseEnhancementFFModel", config)
{

}
//____________________________________________________________________________
TransverseEnhancementFFModel::~TransverseEnhancementFFModel()
{

}
//____________________________________________________________________________
// Return the electric form factor of the base model.
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gep(const Interaction * interaction) const
{
  return fElFormFactorsBase->Gep(interaction);
}
//____________________________________________________________________________
// Return the magnetic form factor of the base model, multiplied by the
// Transverse Enhancement function.
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gmp(const Interaction * interaction) const
{
  return GetTransEnhMagFF(fElFormFactorsBase->Gmp(interaction), interaction);
}
//____________________________________________________________________________
// Return the electric form factor of the base model.
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gen(const Interaction * interaction) const
{
  return fElFormFactorsBase->Gen(interaction);
}
//____________________________________________________________________________
// Return the magnetic form factor of the base model, multiplied by the
// Transverse Enhancement function.
//____________________________________________________________________________
double TransverseEnhancementFFModel::Gmn(const Interaction * interaction) const
{
  return GetTransEnhMagFF(fElFormFactorsBase->Gmn(interaction), interaction);
}
//____________________________________________________________________________
// Multiplies the supplied magnetic form factor by the Transverse Enhancement
// function and returns the result.
//____________________________________________________________________________
double TransverseEnhancementFFModel::GetTransEnhMagFF(
    double magFF, const Interaction * interaction) const
{
  const Target& target = interaction->InitState().Tgt();
  double transEnhA, transEnhB;
  GetTransEnhParams(target, &transEnhA, &transEnhB);
  if (transEnhA == 0) {
    return magFF;
  }
  double Q2 = interaction->Kine().Q2();
  double rt = 1 + transEnhA * Q2 * TMath::Exp(-Q2 / transEnhB);
  return TMath::Sqrt(rt)*magFF;
}
//____________________________________________________________________________
// Returns the Transverse Enhancement parameters as loaded from config files.
//____________________________________________________________________________
void TransverseEnhancementFFModel::GetTransEnhParams(
    const Target& target, double* teA, double* teB) const {
  if (!GetValueFromNuclearMaps(target, fNucMagFF_RT_A,
                               fRangeMagFF_RT_A, teA) ||
      !GetValueFromNuclearMaps(target, fNucMagFF_RT_B,
                               fRangeMagFF_RT_B, teB) ||
      *teB == 0) {
    *teA = 0;
    *teB = 1;
  }
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void TransverseEnhancementFFModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
// Loads Transverse enhancement parameters.  All parameters are from config
// files.
//____________________________________________________________________________
void TransverseEnhancementFFModel::LoadConfig(void)
{
	LoadAllIsotopesForKey("MagFF_RT_A", "TansverseEnhancementFFModel", GetOwnedConfig(),
                        &fNucMagFF_RT_A);
  LoadAllNucARangesForKey("MagFF_RT_A", "TransverseEnhancementFFModel", GetOwnedConfig(),
                          &fRangeMagFF_RT_A);
  LoadAllIsotopesForKey("MagFF_RT_B", "TransverseEnhancementFFModel", GetOwnedConfig(),
                        &fNucMagFF_RT_B);
  LoadAllNucARangesForKey("MagFF_RT_B", "TransverseEnhancementFFModel", GetOwnedConfig(),
                          &fRangeMagFF_RT_B);
}

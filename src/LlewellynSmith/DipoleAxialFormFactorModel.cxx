//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off DipoleELFormFactorsModel by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "LlewellynSmith/DipoleAxialFormFactorModel.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
DipoleAxialFormFactorModel::DipoleAxialFormFactorModel() :
AxialFormFactorModelI("genie::DipoleAxialFormFactorModel")
{

}
//____________________________________________________________________________
DipoleAxialFormFactorModel::DipoleAxialFormFactorModel(string config) :
AxialFormFactorModelI("genie::DipoleAxialFormFactorModel", config)
{

}
//____________________________________________________________________________
DipoleAxialFormFactorModel::~DipoleAxialFormFactorModel()
{

}
//____________________________________________________________________________
double DipoleAxialFormFactorModel::FA(const Interaction * interaction) const
{
  // calculate and return FA
  double q2 = interaction->Kine().q2();
  double fa = fFA0 / TMath::Power(1-q2/fMa2, 2);

  LOG("AxialFormFactor", pDEBUG) << "FA(q^2 = " << q2 << ") = " << fa;
  return fa;
}
//____________________________________________________________________________
void DipoleAxialFormFactorModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleAxialFormFactorModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleAxialFormFactorModel::LoadConfig(void)
{
// get config options from the configuration registry or set defaults 
// from the global parameter list

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // axial mass
  fMa  = fConfig->GetDoubleDef("QEL-Ma", gc->GetDouble("QEL-Ma"));
  fMa2 = TMath::Power(fMa,2);

  // FA(q2 = 0)
  fFA0 = fConfig->GetDoubleDef("QEL-FA0", gc->GetDouble("QEL-FA0"));
}
//____________________________________________________________________________


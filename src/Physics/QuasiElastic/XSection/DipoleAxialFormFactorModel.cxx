//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/DipoleAxialFormFactorModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

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
  GetParam( "QEL-Ma", fMa ) ;
  fMa2 = TMath::Power(fMa,2);

  // FA(q2 = 0)
  GetParam( "QEL-FA0", fFA0 ) ;
}
//____________________________________________________________________________


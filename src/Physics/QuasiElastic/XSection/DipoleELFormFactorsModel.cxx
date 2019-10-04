//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/DipoleELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
DipoleELFormFactorsModel::DipoleELFormFactorsModel() :
ELFormFactorsModelI("genie::DipoleELFormFactorsModel")
{

}
//____________________________________________________________________________
DipoleELFormFactorsModel::DipoleELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::DipoleELFormFactorsModel", config)
{

}
//____________________________________________________________________________
DipoleELFormFactorsModel::~DipoleELFormFactorsModel()
{

}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gep(const Interaction * interaction) const
{
  // calculate and return GNE
  double q2 = interaction->Kine().q2();
  double ge = 1. / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gep(q^2 = " << q2 << ") = " << ge;
  return ge;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gen(const Interaction * /*in*/) const
{
  return 0.;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  // calculate & return Gm
  double q2 = interaction->Kine().q2();
  double gm = fMuP / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmp(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  // calculate & return Gm
  double q2 = interaction->Kine().q2();
  double gm = fMuN / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmn(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::LoadConfig(void)
{
  // vector mass
  GetParam( "EL-Mv", fMv ) ;
  fMv2 = TMath::Power(fMv,2);

  // anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;
}
//____________________________________________________________________________


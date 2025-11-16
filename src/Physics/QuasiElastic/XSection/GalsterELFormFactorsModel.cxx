//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/GalsterELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

//____________________________________________________________________________
GalsterELFormFactorsModel::GalsterELFormFactorsModel() :
ELFormFactorsModelI("genie::GalsterELFormFactorsModel")
{

}
//____________________________________________________________________________
GalsterELFormFactorsModel::GalsterELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::GalsterELFormFactorsModel", config)
{

}
//____________________________________________________________________________
GalsterELFormFactorsModel::~GalsterELFormFactorsModel()
{

}
//____________________________________________________________________________
double GalsterELFormFactorsModel::Gep(const Interaction * interaction) const
{
  double q2  = interaction->Kine().q2();
  double mv2 = fMv2;                         // elastic vector mass^2
  double GD  = 1./TMath::Power(1-q2/mv2,2.); // dipole form factor
  double gep = GD;
  return gep;
}
//____________________________________________________________________________
double GalsterELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  double gmp = fMuP*this->Gep(interaction);
  return gmp;
}
//____________________________________________________________________________
double GalsterELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double q2  = interaction->Kine().q2();
  double M;
  if (fIsIsoscalarNucleon)
  {
    PDGLibrary * pdglib = PDGLibrary::Instance();
    M = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;
  }
  else
  {
      const Target & tgt = interaction->InitState().Tgt();
      M  = tgt.HitNucMass();             // Mnucl
  }

  double M2  = TMath::Power(M,2);            // Mnucl^2
  double t   = -q2/(4*M2);                   // q2<0
  double p   = fGenp;                        // parameter p
  double gen = -1.*fMuN*t*this->Gep(interaction) / (1 + p*t);
  return gen;
}
//____________________________________________________________________________
double GalsterELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  double gmn = fMuN*this->Gep(interaction);
  return gmn;
}
//____________________________________________________________________________
void GalsterELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GalsterELFormFactorsModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void GalsterELFormFactorsModel::LoadConfig(void)
{
  
  GetParamDef( "isIsoscalarNucleon", fIsIsoscalarNucleon, false);
  
  //-- load Galster model parameters

  // Krutov parameters
  GetParam( "Galster-Gen-p", fGenp ) ;

  // vector mass
  GetParam( "EL-Mv",fMv ) ;
  fMv2 = TMath::Power(fMv,2);

  // anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;
}
//____________________________________________________________________________



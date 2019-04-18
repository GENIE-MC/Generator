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
   Renamed LlewellynSmithModel -> LwlynSmithFF
 @ Aug 27, 2013 - AM
   Implemented Axial Form Factor Model structure

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/TransverseEnhancementFFModel.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithFF.h"
#include "Physics/QuasiElastic/XSection/AxialFormFactor.h"
#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF() :
QELFormFactorsModelI()
{

}
//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF(string name) :
QELFormFactorsModelI(name)
{

}
//____________________________________________________________________________
LwlynSmithFF::LwlynSmithFF(string name, string config) :
QELFormFactorsModelI(name, config)
{

}
//____________________________________________________________________________
LwlynSmithFF::~LwlynSmithFF()
{
  if (fCleanUpfElFFModel) {
    delete fElFFModel;
  }
}
//____________________________________________________________________________
double LwlynSmithFF::StrangeF1V(const Interaction * interaction) const
{
  double f1p = this->F1P(interaction); 
  double f1n = this->F1N(interaction);
  double value = 0.;

  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();

  if (pdgc == kPdgSigmaM)        value = -1.* (f1p + 2 * f1n);
  else if (pdgc == kPdgLambda)   value = -kSqrt3 / kSqrt2 * f1p;
  else if (pdgc == kPdgSigma0)   value = -1.* kSqrt2 / 2 * (f1p + 2 * f1n);

  return value;
}
//____________________________________________________________________________
double LwlynSmithFF::StrangexiF2V(const Interaction * interaction) const
{
  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();
  
  double f2p = this->F2P(interaction);
  double f2n = this->F2N(interaction);
  double value = 0.;

  if (pdgc == kPdgSigmaM)
    value = -1.*(f2p +  2.* f2n) ;
  else if (pdgc == kPdgLambda)
    value = (-kSqrt3 / kSqrt2 * f2p) ;
  else if (pdgc == kPdgSigma0)
    value = -1.* kSqrt2 / 2 * (f2p + 2.* f2n) ;

  return value;
}

//____________________________________________________________________________
double LwlynSmithFF::StrangeFA(const Interaction * interaction) const
{
  double value = 0.;

  const XclsTag &      xcls       = interaction->ExclTag();
  int pdgc = xcls.StrangeHadronPdg();

  if (pdgc == kPdgSigmaM)       value =  +1 * (1 - 2 * fFDratio);
  else if (pdgc == kPdgLambda)  value =  -1 / kSqrt6 * (1 + 2 * fFDratio);
  else if (pdgc == kPdgSigma0)  value =  +1 * kSqrt2 / 2 * (1 - 2 * fFDratio);

  fAxFF.Calculate(interaction);
  value *= fAxFF.FA();

  return value;
}
//____________________________________________________________________________
double LwlynSmithFF::F1P(const Interaction * interaction) const
{ 
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gep() - t * fELFF.Gmp());
}
//____________________________________________________________________________
double LwlynSmithFF::F2P(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmp() - fELFF.Gep());
}
//____________________________________________________________________________
double LwlynSmithFF::F1N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gen() - t * fELFF.Gmn());
}
//____________________________________________________________________________
double LwlynSmithFF::F2N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmn() - fELFF.Gen());
}
//____________________________________________________________________________
double LwlynSmithFF::F1V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double _F1V = (gve - t*gvm) / (1-t);
  return _F1V;
}
//____________________________________________________________________________
double LwlynSmithFF::xiF2V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double _xiF2V = (gvm-gve) / (1-t);
  return _xiF2V;
}
//____________________________________________________________________________
double LwlynSmithFF::FA(const Interaction * interaction) const
{
  //-- compute FA(q2) 

  fAxFF.Calculate(interaction);
  return fAxFF.FA();
}
//____________________________________________________________________________
double LwlynSmithFF::Fp(const Interaction * interaction) const
{
  // get momentum transfer
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  // get struck nucleon mass & set pion mass
  const InitialState & init_state = interaction->InitState();
  double MN   = init_state.Tgt().HitNucMass();
  double MN2  = TMath::Power(MN, 2);
  double Mpi  = kPionMass;
  double Mpi2 = TMath::Power(Mpi, 2);

  // calculate FA
  double fa = this->FA(interaction);

  // calculate Fp
  double _Fp = 2. * MN2 * fa/(Mpi2-q2);
  return _Fp;
}
//____________________________________________________________________________
void LwlynSmithFF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithFF::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LwlynSmithFF::LoadConfig(void)
{
// Load configuration data from its configuration Registry (or global defaults)
// to private data members
  fElFFModel =
    dynamic_cast<const ELFormFactorsModelI *> (this->SubAlg("ElasticFormFactorsModel"));
  assert(fElFFModel);

  fCleanUpfElFFModel = false;
  bool useElFFTE = false;
  GetParam( "UseElFFTransverseEnhancement", useElFFTE ) ;
  if(  useElFFTE ) {
    const ELFormFactorsModelI* sub_alg = fElFFModel;
    fElFFModel =
      dynamic_cast<const ELFormFactorsModelI *> (this->SubAlg("TransverseEnhancement"));
    dynamic_cast<const TransverseEnhancementFFModel*>(fElFFModel)->SetElFFBaseModel(
        sub_alg);
    fCleanUpfElFFModel = true;
  }

  fELFF.SetModel(fElFFModel);  

  fAxFFModel =
    dynamic_cast<const AxialFormFactorModelI *> (this->SubAlg("AxialFormFactorModel"));

  assert(fAxFFModel);
  fAxFF.SetModel(fAxFFModel);

  // anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;

  // weinberg angle
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin28w  = TMath::Power(TMath::Sin(thw), 2);

  double d,f ;
  GetParam( "SU3-D", d ) ;
  GetParam( "SU3-F", f ) ;
  fFDratio = f/(d+f); 
}
//____________________________________________________________________________
double LwlynSmithFF::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * MNucl^2)

  //-- get kinematics & initial state parameters
  const Kinematics &   kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  double q2     = kinematics.q2();
  double Mnucl  = init_state.Tgt().HitNucMass();
  double Mnucl2 = TMath::Power(Mnucl, 2);

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*Mnucl2);
}
//____________________________________________________________________________
double LwlynSmithFF::GVE(const Interaction * interaction) const
{
  //-- compute GVE using CVC

  fELFF.Calculate(interaction);
  double gve = fELFF.Gep() - fELFF.Gen();
  return gve;
}
//____________________________________________________________________________
double LwlynSmithFF::GVM(const Interaction * interaction) const
{
  //-- compute GVM using CVC

  fELFF.Calculate(interaction);
  double gvm = fELFF.Gmp() - fELFF.Gmn();
  return gvm;
}
//____________________________________________________________________________


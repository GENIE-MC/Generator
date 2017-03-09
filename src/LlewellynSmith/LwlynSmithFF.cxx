//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
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

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "ElFF/ELFormFactors.h"
#include "ElFF/ELFormFactorsModelI.h"
#include "ElFF/TransverseEnhancementFFModel.h"
#include "Conventions/Constants.h"
#include "LlewellynSmith/LwlynSmithFF.h"
#include "LlewellynSmith/AxialFormFactor.h"
#include "LlewellynSmith/AxialFormFactorModelI.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"

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
	
  fElFFModel = 0;

  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * gc = confp->GlobalParameterList();

  // load elastic form factors model
  RgAlg form_factors_model = fConfig->GetAlgDef(
                 "ElasticFormFactorsModel", gc->GetAlg("ElasticFormFactorsModel"));
  fElFFModel =  dynamic_cast<const ELFormFactorsModelI *> (
                                          this->SubAlg("ElasticFormFactorsModel"));
  assert(fElFFModel);
  fCleanUpfElFFModel = false;
  if(gc->GetBoolDef("UseElFFTransverseEnhancement", false)) {
    const ELFormFactorsModelI* sub_alg = fElFFModel;
    RgAlg transverse_enhancement = gc->GetAlg("TransverseEnhancement");
    fElFFModel = dynamic_cast<const ELFormFactorsModelI *> (
        AlgFactory::Instance()->AdoptAlgorithm(
            transverse_enhancement.name, transverse_enhancement.config));
    dynamic_cast<const TransverseEnhancementFFModel*>(fElFFModel)->SetElFFBaseModel(
        sub_alg);
    fCleanUpfElFFModel = true;
  }

  fELFF.SetModel(fElFFModel);  

  RgAlg ax_form_factor_model = fConfig->GetAlgDef(
                 "AxialFormFactorModel"   , gc->GetAlg("AxialFormFactorModel"   ));

  fAxFFModel =  dynamic_cast<const AxialFormFactorModelI *> (
                                          this->SubAlg("AxialFormFactorModel"   ));
  assert(fAxFFModel);
  fAxFF.SetModel(fAxFFModel);

  // anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));

  // weinberg angle
  double thw = fConfig->GetDoubleDef(
                          "WeinbergAngle", gc->GetDouble("WeinbergAngle"));
  fSin28w  = TMath::Power(TMath::Sin(thw), 2);
  double d = fConfig->GetDoubleDef("QE-SU3-D", gc->GetDouble("SU3-D")); // SU(3) parameter D
  double f = fConfig->GetDoubleDef("QE-SU3-F", gc->GetDouble("SU3-F")); // SU(3) parameter F
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

